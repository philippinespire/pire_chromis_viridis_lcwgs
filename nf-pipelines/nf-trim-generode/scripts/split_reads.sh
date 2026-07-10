#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <in.fastq[.gz]> <out.fastq.gz> <trim_length>" >&2
  exit 1
fi

in="$1"
out="$2"
trim_length="$3"

if [[ "$in" =~ \.gz$ ]]; then
  reader="gunzip -c"
else
  reader="cat"
fi

$reader "$in" | awk -v trim="$trim_length" '
# Line 1: Header
NR%4==1 { header=$0; next }

# Line 2: Sequence
NR%4==2 { sequence=$0; next }

# Line 3: Plus placeholder
NR%4==3 { next }

# Line 4: Quality Scores
NR%4==0 {
  quality=$0;
  L=length(sequence);
  
  # Crucial change: Strict > condition ensures reads <= trim_length are safely passed through intact
  if (L > trim) {
    m=int(L/2); if (m<1) m=1;
    
    if (match(header, /^@[^ ]+/)) {
      id = substr(header, RSTART, RLENGTH);
      rest = substr(header, RSTART+RLENGTH);
      
      # Output First Half
      print id":1" rest;
      print substr(sequence, 1, m);
      print "+";
      print substr(quality, 1, m);
      
      # Output Second Half
      print id":2" rest;
      print substr(sequence, m+1);
      print "+";
      print substr(quality, m+1);
    }
  } else {
    # Keeps original shorter reads perfectly safe, maintaining 4-line sync
    print header;
    print sequence;
    print "+";
    print quality;
  }
}' | gzip -c > "$out"