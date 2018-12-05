#!/usr/bin/env bats

@test "DIAMOND v0.9.23" {
  v="$(diamond --version)"
  [[ "$v" =~ "0.9.23" ]]
}


@test "SRA Toolkit v2.8.2" {
  v="$(fastq-dump --version)"
  [[ "$v" =~ "2.8.2" ]]
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
  python /usr/diamond/helpers/parse_blast.py --input /usr/diamond/tests/unique_sequences.blast --out output.json
  output="$(cat output.json)"
  [[ "$output" =~ "\"id\": \"EcoRI\"" ]]
  [[ "$output" =~ "\"total_depth\": 0.8845" ]]
  [[ "$output" =~ "\"total_reads\": 5" ]]
  [[ "$output" =~ "\"unique_reads\": 5" ]]
  [[ "$output" =~ "\"total_coverage\": 0.8845" ]]
  [[ "$output" =~ "\"unique_coverage\": 0.8845" ]]
}

@test "Parse BLAST results (non-unique reads)" {
  python /usr/diamond/helpers/parse_blast.py --input /usr/diamond/tests/non_unique_sequences.blast --out output.json
  output="$(cat output.json)"
  [[ "$output" =~ "\"id\": \"EcoRI\"" ]]
  [[ "$output" =~ "\"total_depth\": 0.8845" ]]
  [[ "$output" =~ "\"total_reads\": 5" ]]
  [[ "$output" =~ "\"unique_reads\": 4" ]]
  [[ "$output" =~ "\"total_coverage\": 0.8845" ]]
  [[ "$output" =~ "\"unique_coverage\": 0.7076" ]]
}

@test "Make sure the run script is in the PATH" {
  h="$(run.py -h)"

  [[ "$h" =~ "Align a set of reads against a reference database with DIAMOND" ]]
}

@test "Make sure the run_blast.py script is in the PATH" {
  h="$(run_blast.py -h)"

  [[ "$h" =~ "Align a set of reads against a reference database with DIAMOND" ]]
}
