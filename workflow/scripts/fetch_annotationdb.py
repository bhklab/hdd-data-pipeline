import argparse
import json
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional
from urllib.parse import urlparse

import requests
import tqdm

BATCH_SIZE = 50
RETRIES = 3
TIMEOUT = 30


def fetch_json(
	session: requests.Session,
	url: str,
	params: Optional[Dict[str, str]] = None,
	retries: int = RETRIES,
	timeout: int = TIMEOUT,
):
	for attempt in range(retries):
		try:
			resp = session.get(url, params=params, timeout=timeout)
			resp.raise_for_status()
			return resp.json()
		except Exception:
			if attempt == retries - 1:
				raise
			time.sleep(2 * (attempt + 1))


def derive_details_url(db_url: str) -> str:
	parsed = urlparse(db_url)
	path = parsed.path
	if path.endswith("/compound/all"):
		path = path.replace("/compound/all", "/compound/many")
	else:
		path = path.rstrip("/") + "/many"
	return parsed._replace(path=path, query="").geturl()


def chunks(items: List[Dict[str, object]], size: int) -> Iterable[List[Dict[str, object]]]:
	for idx in range(0, len(items), size):
		yield items[idx : idx + size]


def main(db_url: str, output_path: str) -> None:
	outpath = Path(output_path)
	outpath.parent.mkdir(parents=True, exist_ok=True)

	session = requests.Session()
	compound_list = fetch_json(session, db_url)
	if not isinstance(compound_list, list):
		raise ValueError("Expected list response from /compound/all")

	details_url = derive_details_url(db_url)
	missing_cids: List[int] = []
	written = 0

	total_batches = (len(compound_list) + BATCH_SIZE - 1) // BATCH_SIZE
	with outpath.open("w", encoding="utf-8") as handle:
		for batch in tqdm.tqdm(chunks(compound_list, BATCH_SIZE), total=total_batches):
			cids = [str(item.get("cid")) for item in batch if item.get("cid") is not None]
			if not cids:
				continue

			params = {
				"compounds": ",".join(cids),
				"format": "json",
				"bioassay": "true",
				"mechanism": "true",
				"toxicity": "true",
			}
			details = fetch_json(session, details_url, params=params)
			if not isinstance(details, list):
				raise ValueError("Expected list response from /compound/many")

			details_by_cid = {item.get("cid"): item for item in details}
			for drug_info in batch:
				cid = drug_info.get("cid")
				drug_details = details_by_cid.get(cid)
				if drug_details is None:
					if cid is not None:
						missing_cids.append(cid)
					continue
				record = {
					"drug_info": drug_info,
					"drug_details": drug_details,
				}
				handle.write(json.dumps(record, ensure_ascii=True) + "\n")
				written += 1

	if missing_cids:
		print(
			f"Warning: {len(missing_cids)} CIDs missing from /compound/many response",
			flush=True,
		)
	print(f"Wrote {written} compound records to {outpath}", flush=True)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(
		prog="fetch_annotationdb",
		description="Batch fetch AnnotationDB compound details into JSONL",
	)
	parser.add_argument("-u", required=True, help="/compound/all endpoint")
	parser.add_argument("-o", required=True, help="Output JSONL path")
	args = parser.parse_args()

	main(db_url=args.u, output_path=args.o)
