#!/usr/bin/env python

"""Compute a comprehensive set of output .gff3 files from an rpsblast analysis."""

__author__ = "Dr. Thomas van Gurp, Ronald de Jongh MSc, Timo Dijkstra BSc."
__version__ = "1.0"
__title__ = "Compute a comprehensive set of output .gff3 files from an rpsblast analysis."
__author_email__ = "r.d.jongh@naktuinbouw.nl"
__date__ = "20-01-2021"

import argparse
import gzip
import os
from pathlib import Path
from typing import List, Tuple, Dict, Any, Union, Optional
from Bio import SeqIO

# raise NotImplementedError("This script is not done yet.")


def parse_args():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(description='Compute a comprehensive set of output .gff3 files from an rpsblast '
                                                 'analysis.')
    parser.add_argument('--rpbsproc_file', type=Path, help='Path to rpbsproc file')
    parser.add_argument('--inputGFF3', type=Path, help='gff3 input containing protein IDs')
    parser.add_argument('--protein_fasta_path', type=Path, help='protein input for truncated names')
    parser.add_argument('--outpath', type=Path, help='Output path for conserved domains, conserved positions etc')
    parser.add_argument('-db', '--database', type=Path, required=False,
                        help='Path to CDD Database Summary file (cddid_all.tbl)',
                        default=Path('/5_workspace/tools/rpsbproc/RpsbProc-x64-linux/data/cddid.tbl'))
    args = parser.parse_args()
    return args


def get_children(domain_dict):
    """
    Features in GFF3 files can be hierarchical, in that they can have parents and children defined in their ninth
    column. Here we define this hierarchical relationship.
    :param domain_dict: Dictionary containing all entry data for the particular hit.
    :return: child attributes for gff entry.
    """
    i = 0
    regions = []
    # if domain_dict['strand'] == '+':
    start = domain_dict['start']
    # else:
    #     start = domain_dict['end']
    #     # prev_pos = domain_dict['end']
    prev_pos = start - 1
    for i in sorted(domain_dict['mapping_positions']):
        if prev_pos != i - 1:
            stop = prev_pos
            # add child to region
            regions.append((start, stop + 1,))
            start = i
        prev_pos = i
    regions.append((start, i + 1,))
    children = []
    for region in regions:
        child = '%(scaffold)s\tNCBI-CD-batch\tCDS\t' % domain_dict
        child += '%s\t%s\t' % region
        child += '.\t%(strand)s\t' % domain_dict
        # TODO: determine phase
        phase = '0'
        child += '%s\t' % phase
        child += 'Parent=%s_%s' % (domain_dict['name'], domain_dict['protein'])
        children.append(child)
    if not regions:
        # make sure we have a CDS feature for single regions
        child = '%(scaffold)s\tNCBI-CD-batch\tCDS\t' % domain_dict
        child += '%(start)s\t%(end)s\t' % domain_dict
        child += '.\t%(strand)s\t' % domain_dict
        # TODO: determine phase
        phase = '0'
        child += '%s\t' % phase
        child += 'Parent=%s_%s' % (domain_dict['name'], domain_dict['protein'])
        children.append(child)
    return '\n'.join(children)


def format_gff3_entry(domain_dict: Dict[str, Any]) -> str:
    """
    Convert domain in gff3 formatted feature.
    :param domain_dict: Dictionary containing all entry data for the particular hit.
    :return: String of the conserved domains in gff3 format.
    """
    # gff_output will be build item for item
    header = '##sequence-region %(scaffold)s %(start)s %(end)s' % domain_dict
    # gff_string = '{scaffold}\tNCBI-CD-batch\tmRNA\t{start}\t{end}\t{E-value}\t{strand}\t0'.format(domain_dict)
    gff_output = []
    # 1. seqid
    gff_output.append(domain_dict['scaffold'])
    # 2. source
    gff_output.append('NCBI-CD-batch')
    # 3. type
    gff_output.append('mRNA')
    # 4. start
    gff_output.append(str(domain_dict['start']))
    # 5. end
    gff_output.append(str(domain_dict['end']))
    # 6. score
    gff_output.append(domain_dict['E-Value'])
    # 7. strand
    gff_output.append(domain_dict['strand'])
    # 8. phase
    gff_output.append('0')
    # 9 attributes [ID,Name,Alias,Parent,Target,Gap]
    domain_dict['name'] = domain_dict['Short name']  # + '_' + domain_dict['protein']
    domain_dict['E-Value'] = '%.2E' % float(domain_dict['E-Value'])
    if ';' in domain_dict['Definition'][:60]:
        domain_dict['Description'] = '%s EV:%s' % (domain_dict['Definition'].split(';')[0], domain_dict['E-Value'])
    else:
        domain_dict['Description'] = '%s EV:%s' % (domain_dict['Short name'], domain_dict['E-Value'])
    attributes = 'ID=%s_%s' % (domain_dict['name'], domain_dict['protein'])
    translate = {'Definition': 'Long Description',
                 'Short name': None,
                 'E-Value': 'Expectation value (E-value)',
                 'From': 'protein_start_match',
                 'Query': None,
                 'To': 'protein_end_match',
                 'mapping_positions': None,
                 'protein': None,
                 'scaffold': None,
                 'strand': None, }
    for k, v in domain_dict.items():
        if k in translate:
            k = translate[k]
        if k is None or v == ' - ':
            continue
        if type(v) == str:
            v = v.replace(';', '%3B').replace(',', '%2C').replace('=', '%3D').replace('&', '%26')
        attributes += ';%s=%s' % (k, v)
    attributes += ';Alias=%s,' % domain_dict['Accession']
    if domain_dict['Superfamily'] != ' - ':
        attributes += '%s,' % domain_dict['Superfamily']
    attributes += '%s,' % domain_dict['PSSM-ID']
    attributes += '%s' % domain_dict['Description']
    gff_output.append(attributes)
    # make child attributes for CDS features
    children = get_children(domain_dict)
    if children:
        output = '\n'.join([header, '\t'.join(gff_output), children]) + '\n'
    else:
        output = header + '\n' + '\t'.join(gff_output) + '\n'
    return output


def format_gff3_conspos(domain_dict: Union[Dict[str, Any], None]) -> str:
    """
    Convert conserved positions in gff3 formatted feature.
    :param domain_dict: Dictionary containing all entry data for the particular hit.
    :return: String of the conserved positions in gff3 format.
    """
    # gff_output will be build item for item
    header = '##sequence-region %(scaffold)s %(start)s %(end)s' % domain_dict
    gff_output = []
    # 1. seqid
    gff_output.append(domain_dict['scaffold'])
    # 2. source
    gff_output.append('NCBI-CD-batch')
    # 3. type
    gff_output.append('mRNA')
    # 4. start
    gff_output.append(str(domain_dict['start']))
    # 5. end
    gff_output.append(str(domain_dict['end']))
    # 6. score
    gff_output.append(domain_dict['mapped size'])
    # 7. strand
    gff_output.append(domain_dict['strand'])
    # 8. phase
    gff_output.append('0')
    # 9 attributes [ID,Name,Alias,Parent,Target,Gap]
    domain_dict['name'] = domain_dict['Title']
    domain_dict['Description'] = domain_dict['Title']
    attributes = 'ID=%s_%s' % (domain_dict['name'], domain_dict['protein'])
    translate = {'Query': None,
                 'Title': None,
                 'mapping_positions': None,
                 'protein': None,
                 'scaffold': None,
                 'strand': None, }
    for k, v in domain_dict.items():
        if k in translate:
            k = translate[k]
        if k is None or v == ' - ':
            continue
        if type(v) == str:
            v = v.replace(';', '%3B').replace(',', '%2C').replace('=', '%3D').replace('&', '%26')
        attributes += ';%s=%s' % (k, v)
    gff_output.append(attributes)
    # make child attributes for CDS features
    children = get_children(domain_dict)
    if children:
        output = '\n'.join([header, '\t'.join(gff_output), children]) + '\n'
    else:
        output = header + '\n' + '\t'.join(gff_output) + '\n'
    return output


def parse_gff3(path: str) -> dict:
    """
    parse gff3 object
    :param path: Path to the .gff3 file as string.
    :return: A dictionary containing positional metadata for every protein: scaffold, strand, start, end, positions.
    """
    mapping_dict = {}
    if path.endswith('.gz'):
        handle = gzip.open(path)
    else:
        handle = open(path)
    line = handle.readline()
    for line in handle:
        if line.strip() == '##FASTA':
            break
        if line.startswith('#'):
            continue
        if line.startswith('A') or line.startswith('T') or line.startswith(
                'C') or line.startswith('G') or line.startswith('>'):
            print('DNA between normal lines!')
            continue
        split_line = line.rstrip('\n').split('\t')
        if split_line[2] in ['mRNA', 'CDS']:
            scaffold, source, type_, start, end, score, strand, phase, attributes = split_line
            if type_ == 'CDS':
                # phase determines start/end position of exon start amino acid in CDS.
                start = int(start)
                end = int(end)
                if strand == '+':
                    if phase == '1':
                        start += 1
                    elif phase == '2':
                        start -= 1
                    # Then end of the CDS does not necesarily correspond to the last amino acid in the exon
                    # take modulo to determine end of amino acid in CDS frame
                    mod = (end - start) % 3
                    if mod == 0:
                        end -= 1
                    elif mod == 1:
                        end += 1
                    pass
                elif strand == '-':
                    if phase == '1':
                        end -= 1
                    elif phase == '2':
                        end += 1
                    # The start of the CDS does not necesarily correspond to the last amino acid in the exon
                    # take modulo to determine start of amino acid in CDS frame
                    mod = (end - start) % 3
                    if mod == 0:
                        start -= 1
                    elif mod == 1:
                        start += 1
                    pass
                else:
                    continue
            if type_ == 'mRNA':
                id = ''
                index_word = 'ID='
                for nt in attributes[attributes.index(index_word) + len(index_word):]:
                    if nt not in [';', '\n', ]:
                        id += nt
                    else:
                        break
                if id.startswith('mRNA:'):
                    id = id[5:]
            else:
                # fixme: this code will never be reached?
                id = ''
                index_word = 'ID='
                for nt in attributes[attributes.index(index_word) + len(index_word):]:
                    if nt not in [';', '\n', ]:
                        id += nt
                    else:
                        break
                if id.startswith('mRNA:'):
                    id = id[5:]
                # id = ''
                # for nt in attributes[attributes.index('Parent=') + len('Parent='):]:
                #     if nt not in [';','\n',]:
                #         id += nt
                #     else:
                #         break
                # if id.startswith('mRNA:'):
                #     id = id[5:]
            if type_ == 'mRNA':
                mapping_dict[id] = {'scaffold': scaffold, 'strand': strand, 'start': int(start),
                                    'end': int(end)}
            else:
                if id not in mapping_dict:
                    mapping_dict[id] = {'scaffold': scaffold, 'strand': strand, 'start': int(start),
                                        'end': int(end)}
                positions = range(int(start), int(end + 1))
                try:
                    mapping_dict[id]['positions'] += positions
                except KeyError:
                    mapping_dict[id]['positions'] = positions
    handle.close()
    # now a protein position can be mapped back to the reference position
    return mapping_dict


def protein_index(protein: Path) -> Dict[int, str]:
    """
    Create a dictionary with protein IDs.
    :param protein: protein input for truncated names.
    :return: Dictionary in the format key: protein index, value: protein ID.
    """
    with open(protein) as handle:
        seqrecords = [record for record in SeqIO.parse(handle, "fasta")]
    nr_id = {index + 1: record.id for index, record in enumerate(seqrecords)}
    return nr_id


def parse_domain(line: str, header: List[str], mapping_dict: dict, nr_id: Dict[int, str]) -> Tuple[str, Dict[str, Any]]:
    """
    Parse domain into object
    :param line: The line of interest from hitdata.txt.
    :param header: List of the column names from hitdata.txt.
    :param mapping_dict: A dictionary containing positional metadata for every protein: scaffold, strand, start, end,
    positions.
    :param nr_id: Dictionary in the format key: protein index, value: protein ID.
    :return: The entry as string formatted in GFF3 format: the entry header, parent, child. In addition, a dictionary is
    returned containing all entry data.
    """
    domain_dict = {}
    for h, i in zip(header, line[:-1].split('\t')):
        if h in ['domainStarts', 'domainEnds']:
            i = i.split(',')[0]
        domain_dict[h] = i
    nr = int(domain_dict['Query'][2:domain_dict['Query'].index(' ')])
    try:
        domain_dict['protein'] = nr_id[nr].lstrip('>')
    except KeyError:
        pass
    id = domain_dict['protein']
    if id in mapping_dict:
        prot_start = int(domain_dict['From']) * 3
        prot_end = int(domain_dict['To']) * 3

        if prot_end >= len(mapping_dict[id]['positions']):
            prot_end = len(mapping_dict[id]['positions']) - 1
        if mapping_dict[id]['strand'] == '+':
            positions = sorted(mapping_dict[id]['positions'])
        else:
            positions = sorted(mapping_dict[id]['positions'])[::-1]
        seq_start, seq_end = [positions[prot_start], positions[prot_end]]
        # determine seq-range
        mapping_positions = []
        if seq_start < seq_end:
            pos_range = range(seq_start, seq_end)
            if not pos_range:
                pos_range = range(seq_end, seq_start)
        elif seq_start > seq_end:
            pos_range = range(seq_end, seq_start)
            if not pos_range:
                pos_range = range(seq_start, seq_end)
        # print(mapping_dict[id]['positions'])
        # print(pos_range)
        # print('\n')
        for pos in pos_range:
            if pos in mapping_dict[id]['positions']:
                mapping_positions.append(pos)
        # print(mapping_positions)
        domain_dict['mapping_positions'] = mapping_positions
        # print(mapping_positions[0])
        domain_dict['start'] = mapping_positions[0]
        domain_dict['end'] = mapping_positions[-1]
        scaffold = mapping_dict[id]['scaffold']
        domain_dict['scaffold'] = scaffold
        if seq_start > seq_end:
            strand = '-'
        else:
            strand = '+'
        domain_dict['strand'] = strand
        entry = format_gff3_entry(domain_dict)
        return entry, domain_dict
    else:
        raise KeyError('id %s not found in mapping dict' % id)


def parse_conserved_positions(line: str, header: List[str], mapping_dict: dict, nr_id: Dict[int, str]) -> \
        Union[Tuple[str, Dict[str, Any]], None]:
    """
    Parse domain into object
    :param line: The line of interest from featdata.txt.
    :param header: List of the column names from featdata.txt.
    :param mapping_dict: A dictionary containing positional metadata for every protein: scaffold, strand, start, end,
    positions.
    :param nr_id: Dictionary in the format key: protein index, value: protein ID.
    :return: The entry as string formatted such that it starts with a header and lists the conserved site on lines
    below in GFF3 format. In addition, a dictionary is returned containing all entry data.
    """
    domain_dict = {}
    for h, i in zip(header, line[:-1].split('\t')):
        if h in ['domainStarts', 'domainEnds']:
            i = i.split(',')[0]
        domain_dict[h] = i
    nr = int(domain_dict['Query'][2:domain_dict['Query'].index(' ')])
    try:
        domain_dict['protein'] = nr_id[nr].lstrip('>')
    except KeyError:
        pass
    id = domain_dict['protein']
    if id in mapping_dict:
        # prot_positions contains the locations of conserved features on the protein
        prot_positions = []
        if mapping_dict[id]['strand'] == '+':
            positions = sorted(mapping_dict[id]['positions'])
            offset = 1
        else:
            positions = sorted(mapping_dict[id]['positions'])[::-1]
            offset = -1
        if ',' in domain_dict['coordinates']:
            for loc in domain_dict['coordinates'].split(','):
                aa = loc[0]
                pos = int(loc[1:]) - 1
                prot_positions.append((3 * pos))
            mapping_positions = []
            for prot_pos in prot_positions:
                try:
                    if mapping_dict[id]['strand'] == '+':
                        nt_pos = positions[prot_pos]
                    else:
                        nt_pos = positions[prot_pos + 1]
                    mapping_positions += [nt_pos, nt_pos + offset]
                except IndexError:
                    break
        elif '-' in domain_dict['coordinates']:
            # range is given for converved nucleotides in repeat domain
            start, end = domain_dict['coordinates'].split('-')
            start_pos = 3 * (int(start[1:]) - 1)
            end_pos = 3 * (int(end[1:]) - 1)
            # we need sort because otherwise range won't work
            mapping_positions = sorted([positions[start_pos], positions[end_pos + (offset * 2)]])
            mapping_positions = range(mapping_positions[0], mapping_positions[1])
        else:
            # single position
            loc = domain_dict['coordinates']
            aa = loc[0]
            pos = 3 * (int(loc[1:]) - 1)
            mapping_positions = sorted([positions[pos], positions[pos + (offset * 2)]])
            mapping_positions = range(mapping_positions[0], mapping_positions[1])
        scaffold = mapping_dict[id]['scaffold']
        domain_dict['scaffold'] = scaffold
        domain_dict['strand'] = mapping_dict[id]['strand']
        domain_dict['mapping_positions'] = mapping_positions
        domain_dict['start'] = min(mapping_positions[0], mapping_positions[-1])
        domain_dict['end'] = max(mapping_positions[0], mapping_positions[-1]) + 1
        entry = format_gff3_conspos(domain_dict)
        return entry, domain_dict
    else:
        return None


def conserved_domains(hitdata_path: Path, protein_fasta_path: Path, gff_input: Path, outpath: Path):
    """
    Go over the hitdata and start building .gff3 output files by reading a line from hitdata, converting it to .gff3
    and appending it to the corresponding output file.
    :param hitdata_path: Path to the hitdata.txt file.
    :param protein_fasta_path: protein input for truncated names.
    :gff_input: gff3 input containing protein IDs.
    :param outpath: Output path for conserved domains, conserved positions etc.
    """
    if hitdata_path.suffix == '.gz':
        hitdata = gzip.open(str(hitdata_path))
    else:
        hitdata = hitdata_path.open()
    nr_id = protein_index(protein_fasta_path)

    while True:
        header = hitdata.readline()
        if header.startswith('Query'):
            header = header[:-1].split('\t')
            break

    gff_header = "##gff-version 3.2.1\n"

    mapping_dict = parse_gff3(str(gff_input))
    for line in hitdata:
        entry, domain_dict = parse_domain(line, header, mapping_dict, nr_id)
        out_path = os.path.join('%s' % outpath,
                                '%s.gff3' % domain_dict['Hit type'])
        with open(out_path, 'a') as file_out:
            file_out.write(entry)


def conserved_sites(featdata: Path, protein_fasta_path: Path, gff_input: Path, outpath: Path):
    """
    Go over the featdat and start building .gff3 output files by reading a line from featdata, converting it to .gff3
    and appending it to the corresponding output file.
    :param featdata: Path to the featdata.txt file.
    :param protein_fasta_path: protein input for truncated names.
    :gff_input: gff3 input containing protein IDs.
    :param outpath: Output path for conserved domains, conserved positions etc.
    """
    out_path = os.path.join(outpath, 'conserved_sites.gff3')
    with open(out_path, 'a') as file_out:

        if featdata.suffix == 'gz':
            featdata = gzip.open(featdata)
        else:
            featdata = featdata.open()
        nr_id = protein_index(protein_fasta_path)
        while True:
            header = featdata.readline()
            if header.startswith('Query'):
                header = header[:-1].split('\t')
                break
        mapping_dict = parse_gff3(str(gff_input))
        for line in featdata:
            try:
                entry, domain_dict = parse_conserved_positions(line, header, mapping_dict, nr_id)
            except TypeError:
                continue
            file_out.write(entry)


def process_data(list_of_lines: List[str], startmotif: str, endmotif: str, db: Path, trigger: bool) -> List[str]:
    """
    Read the .rpsbproc file. Filter input file by motif, either by domain or site. If the startmotif matches the desired
    input, retrieve results. In addition, get a domain description from the database if the trigger is true.
    :param list_of_lines: List of lines from the .rpsbproc file.
    :param startmotif: Either <SITES> or <DOMAINS>.
    :param endmotif: Either <ENDSITES> or <ENDDOMAINS>.
    :param db: Path to the cddid.tbl file.
    :param trigger: Boolean, should be True if you process the data for domains.
    :return: List of results from the .rpsbproc file.
    """

    data = []
    copy = False
    query_definition = ''
    for line in list_of_lines:
        if not line.startswith('#'):
            if line.strip().startswith('QUERY'):
                query_definition = line.split('\t')[4]
            elif line.strip() == startmotif:
                copy = True
            elif line.strip() == endmotif:
                copy = False
            elif copy:
                line = line.strip('\n').split('\t')[1:]
                line[0] = line[0][0] + '#' + line[0][6:] + ' - >' + query_definition.strip()
                if trigger:
                    line.append(domain_description(line[2], db))
                data.append('\t'.join(line))
    return data


# def domain_description(PSSMid):
#     url = 'https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid={}'.format(PSSMid)
#     page = requests.get(url, allow_redirects=True)
#     soup = BeautifulSoup(page.text, 'html.parser')
#     description = soup.find('div', attrs={'class': 'inner'})
#     try:
#         description = description.text.strip()
#         if '\n' in description:
#             description = description.split('\n')[4]
#     except AttributeError:
#         description = '-'
#     print(description)
#     return description


def domain_description(PSSMid: str, db: Path) -> str:
    """
    Search for domain descriptions in local CDD database with PSSMid.
    :param PSSMid: PSSMid as string.
    :param db: Path to the cddid.tbl file.
    :return: Domain description corresponding to the PSSMid.
    """

    with open(db, ) as d:
        line = next((l for l in d if PSSMid in l), None)
    return line.split('\t')[3]


def processes_rpsbproc_output(input_path: Path, database: Path, outpath) -> Tuple[Path, Path]:
    """
    Read the .rpsbproc file and generate featdata for the sites and hitdata for the domains. Write this to two
    output .txt files.
    :param input_path: Path to rpbsproc file.
    :param database: Path to CDD Database Summary file (cddid_all.tbl).
    :param outpath: Output path for conserved domains, conserved positions etc.
    :return: Tuple of the paths to the featdata and hitdata text files.
    """
    feat_header = 'Query\tType\tTitle\tcoordinates\tcomplete size\tmapped size\tsource domain\n'
    hit_header = 'Query\tHit type\tPSSM-ID\tFrom\tTo\tE-Value\tBitscore\tAccession\tShort name\tIncomplete\t' \
                 'Superfamily\tDefinition\n'

    with input_path.open() as f:
        featdata = process_data(f.readlines(), 'SITES', 'ENDSITES', database, False)
        f.seek(0)
        hitdata = process_data(f.readlines(), 'DOMAINS', 'ENDDOMAINS', database, True)

    featdata_path = (outpath / 'featdata.txt')
    featdata_path.write_text(feat_header + '\n'.join(featdata))
    hitdata_path = (outpath / 'hitdata.txt')
    hitdata_path.write_text(hit_header + '\n'.join(hitdata))
    return featdata_path, hitdata_path


def main(rpbsproc_file: Path, outpath: Path, protein_fasta_path: Path, inputGFF3: Path,
         database: Optional[Path] = Path('/5_workspace/tools/rpsbproc/RpsbProc-x64-linux/data/cddid.tbl')):
    featdata_path, hitdata_path = processes_rpsbproc_output(rpbsproc_file, database, outpath)

    conserved_domains(hitdata_path, protein_fasta_path, inputGFF3, outpath)
    conserved_sites(featdata_path, protein_fasta_path, inputGFF3, outpath)


if __name__ == '__main__':
    args = parse_args()
    ext_list = "\n".join([str(_) for _ in (args.rpbsproc_file, args.protein_fasta_path, args.inputGFF3,
                                      args.database) if not _.exists()])
    assert not ext_list, f'Did not exist:\n{ext_list}'
    if args.outpath.exists():
        print(f'Warning, file exists: {args.outpath}')
    main(args.rpbsproc_file, args.outpath, args.protein_fasta_path, args.inputGFF3, database=args.database)
