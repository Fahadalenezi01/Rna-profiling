# filter_gtf_by_id.py
import sys
import re
import argparse

def main():
    parser = argparse.ArgumentParser(description="Filter GTF file for transcript lines matching IDs in a list.")
    parser.add_argument("-i", "--ids", required=True, help="File containing transcript IDs (one per line).")
    parser.add_argument("-g", "--gtf", required=True, help="Input GTF file to filter.")
    parser.add_argument("-o", "--output", required=True, help="Output GTF file.")
    args = parser.parse_args()

    print(f"Reading IDs from: {args.ids}")
    try:
        with open(args.ids, 'r') as f_ids:
            # Read IDs and store them in a set for fast lookup
            target_ids = set(line.strip() for line in f_ids if line.strip())
        print(f"Read {len(target_ids)} unique target IDs.")
        if not target_ids:
            print("Warning: ID list file is empty.")
            # Create empty output file and exit
            open(args.output, 'w').close()
            sys.exit(0)
    except FileNotFoundError:
        print(f"Error: ID file not found: {args.ids}")
        sys.exit(1)

    print(f"Filtering GTF file: {args.gtf}")
    written_count = 0
    try:
        with open(args.gtf, 'r') as f_gtf, open(args.output, 'w') as f_out:
            for line in f_gtf:
                # Skip comments
                if line.startswith('#'):
                    # Optionally write comments to output: f_out.write(line)
                    continue

                fields = line.strip().split('\t')
                # Ensure it's a valid GTF line with at least 9 fields
                if len(fields) < 9:
                    continue

                # Check if it's a transcript line
                if fields[2].lower() == 'transcript':
                    attributes = fields[8]
                    # Use regex to find transcript_id "VALUE";
                    match = re.search(r'transcript_id "([^"]+)"', attributes)
                    if match:
                        transcript_id = match.group(1)
                        # Check if this ID is in our target set
                        if transcript_id in target_ids:
                            f_out.write(line)
                            written_count += 1
    except FileNotFoundError:
        print(f"Error: Input GTF file not found: {args.gtf}")
        sys.exit(1)
    except Exception as e:
         print(f"An error occurred during GTF processing: {e}")
         sys.exit(1)

    print(f"Filtering complete. Wrote {written_count} matching transcript lines to: {args.output}")
    if written_count == 0:
         print("Warning: No matching transcript lines found in the input GTF.")

if __name__ == "__main__":
    main()
