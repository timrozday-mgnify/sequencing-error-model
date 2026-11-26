#!/usr/bin/env python3
"""Fetch ENA files from FIRE's S3 endpoint into the working directory, under their ENA names.

As raw-reads-analysis-pipeline's s3fire_downloader.py: ENA's vol1/ is the era-public bucket (read unsigned) or, for
private data, era-private (signed with the FIRE_ACCESS_KEY / FIRE_SECRET_KEY Nextflow secrets).
"""

import os
import re
import sys

import boto3
from botocore import UNSIGNED
from botocore.config import Config

ENDPOINT = "https://hl.fire.sdo.ebi.ac.uk"
BUCKETS = {"ftp.sra.ebi.ac.uk": "era-public", "ftp.dcc-private.ebi.ac.uk": "era-private"}


def fire_key(url: str) -> tuple[str, str]:
    """The bucket and key of an ENA vol1 URL (ftp://, http:// or https://, or no scheme)."""
    m = re.match(r"^(?:\w+://)?([^/]+)/vol1/(.+)$", url)
    if not m or m[1] not in BUCKETS:
        raise ValueError(f"not an ENA vol1 URL: {url}")
    return BUCKETS[m[1]], m[2]


def main(urls: list[str]) -> None:
    for url in urls:
        bucket, key = fire_key(url)
        if bucket == "era-public":
            s3 = boto3.client("s3", endpoint_url=ENDPOINT, config=Config(signature_version=UNSIGNED))
        else:
            s3 = boto3.client(
                "s3",
                endpoint_url=ENDPOINT,
                aws_access_key_id=os.environ["FIRE_ACCESS_KEY"],
                aws_secret_access_key=os.environ["FIRE_SECRET_KEY"],
            )
        print(f"{url} -> s3://{bucket}/{key}", file=sys.stderr)
        s3.download_file(bucket, key, os.path.basename(key))


if __name__ == "__main__":
    main(sys.argv[1:])
