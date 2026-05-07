#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <in.fastq[.gz]> <out.fastq.gz> <trim_length>" >&2
  exit 1
fi

in="$1"
out="$2"
trim_length="$3"

reader="cat"
[[ "$in" =~ \.gz$ ]] && reader="zcat"

$reader "$in" | awk 'NR%4==1{
  header=$0; getline sequence; getline plus; getline quality;
  L=length(sequence);
  if (L>="'"$trim_length"'") {
    m=int(L/2); if (m<1) m=1;
    # Split header: @ID rest_of_header
    # Append :1 and :2 to the ID portion
    if (match(header, /^@[^ ]+/)) {
      id = substr(header, RSTART, RLENGTH);
      rest = substr(header, RSTART+RLENGTH);
      # first half
      print id":1" rest;
      print substr(sequence,1,m);
      print plus;
      print substr(quality,1,m);
      # second half
      print id":2" rest;
      print substr(sequence,m+1);
      print plus;
      print substr(quality,m+1);
    }
  } else {
    print header; print sequence; print plus; print quality;
  }
}' | gzip -c > "$out"