#!/usr/bin/env bats

@test "DIAMOND v0.9.10" {
  v="$(diamond --version)"
  [[ "$v" =~ "0.9.10" ]]
}

@test "Parse BLAST results" {
  python /usr/diamond/helpers/parse_blast.py --input /usr/diamond/tests/reads.blast --out output.json
  output="$(cat output.json)"
  [[ "$output" =~ "\"id\": \"EcoRI\"" ]]
  [[ "$output" =~ "\"total_depth\": 0.8845" ]]
  [[ "$output" =~ "\"total_coverage\": 0.8845" ]]
}

@test "Make sure the run script is in the PATH" {
  h="$(run.py -h)"

  [[ "$h" =~ "Align a set of reads against a reference database with DIAMOND" ]]
}
