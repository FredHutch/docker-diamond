#!/usr/bin/env bats

@test "DIAMOND v0.9.10" {
  v="$(diamond --version)"
  [[ "$v" =~ "0.9.10" ]]
}

@test "AWS CLI v1.11.146" {
  v="$(aws --version 2>&1)"
  [[ "$v" =~ "1.11.146" ]]
}


@test "Curl v7.47.0" {
  v="$(curl --version)"
  [[ "$v" =~ "7.47.0" ]]
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
