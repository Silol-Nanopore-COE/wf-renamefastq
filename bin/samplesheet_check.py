#!/usr/bin/env python

# This script is created by the combination of functions in the sample sheet check script from
# epi2me-labs: wf-bacterial-genomes -> https://github.com/epi2me-labs/wf-bacterial-genomes/blob/master/bin/workflow_glue/check_sample_sheet.py
# nf-core: viralrecon -> https://github.com/nf-core/viralrecon/blob/master/bin/check_samplesheet.py

import codecs
import sys
import re
import argparse


def parse_args(args=None):
    Description = "Check if a sample sheet is valid."
    Epilog = "Example usage: python samplesheet_check.py <SAMPLE_SHEET>"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-i",
        "--sample_sheet",
        dest="SAMPLE_SHEET",
        help="Input samplesheet file in the CSV file format.",
    )
    return parser.parse_args(args)


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


# Some Excel users save their CSV as UTF-8 (and occasionally for a reason beyond my
# comprehension, UTF-16); Excel then adds a byte order mark (unnecessarily for UTF-8
# I should add). If we do not handle this with the correct encoding, the mark will
# appear in the parsed data, causing the header to be malformed.
# See CW-2310
def determine_codec(file_in):
    """Peek at a file and return an appropriate reading codec."""
    with open(file_in, "rb") as f_bytes:
        # Could use chardet here if we need to expand codec support
        initial_bytes = f_bytes.read(8)

        for codec, encoding_name in [
            [codecs.BOM_UTF8, "utf-8-sig"],  # use the -sig codec to drop the mark
            [codecs.BOM_UTF16_BE, "utf-16"],  # don't specify LE or BE to drop mark
            [codecs.BOM_UTF16_LE, "utf-16"],
            [codecs.BOM_UTF32_BE, "utf-32"],  # handle 32 for completeness
            [codecs.BOM_UTF32_LE, "utf-32"],  # again skip LE or BE to drop mark
        ]:
            if initial_bytes.startswith(codec):
                return encoding_name
        return None  # will cause file to be opened with default encoding


def samplesheet_check(sample_sheet):
    """
    This function checks that the samplesheet follows the following structure:

    barcode,alias
    barcode01,sample1
    barcode02,sample2

    """
    encoding = determine_codec(sample_sheet)
    with open(sample_sheet, "r", encoding=encoding) as f:
        ## Check header
        MIN_COLS = 2
        HEADER = ["barcode", "alias"]
        header = [x.strip('"') for x in f.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(
                "ERROR: Please check samplesheet header -> {} != {}".format(
                    ",".join(header), ",".join(HEADER)
                )
            )
            sys.exit(1)

        ## Check sample entries
        barcodes = []
        aliases = []
        for line in f:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(
                        MIN_COLS
                    ),
                    "Line",
                    line,
                )
            ## Check barcode and alias entry
            barcode, alias = lspl[: len(HEADER)]

            # Check barcode entry
            if barcode:
                if not re.match(r"^barcode\d\d+$", barcode):
                    print_error(
                        "Barcode entry is not in the correct format!", "Line", line
                    )
            barcodes.append(barcode)

            # Check alias entry
            if alias:
                if alias.find(" ") != -1:
                    print_error(
                        "Spaces have been detected in alias entry: {}".format(alias),
                        "Line",
                        line,
                    )
            aliases.append(alias)

        ## Check if barcode entry is unique in the barcode column
        if len(barcodes) != len(set(barcodes)):
            print_error(
                "Samplesheet contains duplicate entries in the 'barcode' column!",
                "Line",
                line,
            )
        if len(aliases) != len(set(aliases)):
            print_error(
                "Samplesheet contains duplicate entries in the 'alias' column!",
                "Line",
                line,
            )


def main(args=None):
    args = parse_args(args)
    samplesheet_check(args.SAMPLE_SHEET)


if __name__ == "__main__":
    sys.exit(main())
