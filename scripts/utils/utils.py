import argparse
from colorama import Fore
from typing import TextIO, Union, List
import os
from pathlib import Path
import pandas as pd
import bs4.element
from bs4 import BeautifulSoup
import cssutils

NO_NCBI_GENOME_ENTRY = "NO NCBI ENTRY"
MAX_CHARACTERS_HTML_LINE = 9900
IDS_PER_LINE_AFTER_SPLIT = 100


def is_valid_file(parser: argparse, arg: str) -> TextIO:
    """
    Check if the input file exists.
    :param parser: Command line arguments.
    :param arg: The input path.
    :return: If the input path exists, return the opened file.
    """
    if not os.path.isfile(arg):
        parser.error(F"{Fore.RED}The file {arg} does not exist")
    else:
        return open(arg, "r")


def is_valid_directory(parser: argparse, arg: str) -> str:
    """
    Check if the output directory exists.
    :param parser: Command line arguments.
    :param arg: The output path.
    :return: If the input path exists, return the given argument.
    """
    if not os.path.isdir(arg):
        parser.error(F"{Fore.RED}The directory {arg} does not exist")
    else:
        return arg


def is_valid_path(parser: argparse, arg: str) -> str:
    """
    Check if the input file exists.
    :param parser: Command line arguments.
    :param arg: The input path.
    :return: If the input path exists, return the given argument.
    """
    if not os.path.isfile(arg):
        parser.error(F"{Fore.RED}The file {arg} does not exist")
    else:
        return arg


def directory_exists(directory: Path):
    """
    Check if the directory exists. If not, quit.
    :param directory: The output path.
    """
    if not directory.is_dir():
        print(F"{Fore.RED}The directory {directory} does not exist")
        quit()


def file_exists(file: Path):
    """
    Check if the file exists. If not, quit.
    :param file: The path to the file.
    """
    if not file.is_file():
        print(F"{Fore.RED}The file {file} does not exist")
        quit()


def any_effectors(t3sepp_df: pd.DataFrame) -> bool:
    """
    Quit the program if no effector is found by T3SEpp.
    :param t3sepp_df: Dataframe converted T3SEpp output.
    :return: True if no effectors are present.
    """
    return (t3sepp_df.values == "T3S").sum() == 0


def make_dataframe_effectors(t3sepp_out: Union[TextIO, Path]) -> pd.DataFrame:
    """
    Function that takes only the proteins from the T3SEpp output that are predicted as effectors.
    :param t3sepp_out: Direct output from T3SEpp (T3SEpp.out.txt file).
    :return: Dataframe containing only proteins predicted as effector.
    """
    t3sepp_df = pd.read_csv(t3sepp_out, sep="\t")
    if any_effectors(t3sepp_df):
        print(F"{Fore.BLUE}No effectors found.")
        quit()

    t3sepp_df_effectors = t3sepp_df.loc[t3sepp_df["Pred"] == "T3S"]
    return t3sepp_df_effectors


def read_css(style_tag: bs4.element.Tag) -> dict:
    """
    Scan the html file using BeautifulSoup. Retrieve the styling tag. Convert the BeautifulSoup object to byte format
    and read this by cssutils. Save the css selectors in a dictionary.
    :param style_tag: style tag of the BeautifulSoup object.
    :return: Dictionary with as key the selector and as value the rule.
    """
    selectors = {}
    css = cssutils.parseString(style_tag.encode_contents())
    for rule in css:
        if rule.type == rule.STYLE_RULE:
            style = rule.selectorText
            selectors[style] = {}
            for item in rule.style:
                property_name = item.name
                value = item.value
                selectors[style][property_name] = value
    return selectors


def split_css(selectors: dict) -> List[dict]:
    """
    Determine the length of the css selector. If the length exceeds 9900, split the rule into multiple rules.
    :param selectors: Dictionary with as key the selector and as value the rule.
    :return: List of selector dictionaries.
    """
    selectors_split = []
    for id, rule in selectors.items():
        selector_split = {}
        if len(id) > MAX_CHARACTERS_HTML_LINE:
            id_split = id.split(",")
            for i in range(0, len(id_split), IDS_PER_LINE_AFTER_SPLIT):
                selector_split[", ".join(id_split[i:i+IDS_PER_LINE_AFTER_SPLIT])] = rule
        else:
            selector_split[id] = rule
        selectors_split.append(selector_split)
    return selectors_split


def deal_with_large_tables(output_file: Path) -> BeautifulSoup:
    """
    Pandas has a bug in which css styling has a maximum of 10,000 characters per line. So large dataframes will be cut-
    off at a certain point with respect to styling. Therefore, we came up with a work around that reads the html file
    and will refill the style tag with less ids per line.
    :param output_file: Path to .html output file.
    :return BeautifulSoup object.
    """
    with open(output_file) as html:
        html = html.read()
        soup = BeautifulSoup(html, "html.parser")
    style_tag = soup.select("style")[0]
    selectors = read_css(style_tag)
    selectors_split = split_css(selectors)
    style_tag.clear()
    for selector in selectors_split:
        for key, rule in selector.items():
            rule = str(rule).replace("'", "").replace("{", "{\n").replace(",", ";\n").replace("}", ";}") + "\n"
            css = "{id}{rule}".format(id=key, rule=rule)
            style_tag.append(css)
    return soup
