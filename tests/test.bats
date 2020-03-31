#!/usr/bin/env bats

@test "DIAMOND v0.9.31" {
  v="$(diamond --version)"
  [[ "$v" =~ "0.9.31" ]]
}
