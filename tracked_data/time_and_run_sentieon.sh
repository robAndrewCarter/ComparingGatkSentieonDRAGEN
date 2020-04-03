#!/usr/bin/env bash
JSON_FILE=$1
LOG_FILE=${JSON_FILE}.run.log
echo Start: $(date +'%m-%d-%Y-%H-%M-%S') > $LOG_FILE
echo Running: ~/sentieon-google-genomics/runner/sentieon_runner.py $JSON_FILE
echo $(~/sentieon-google-genomics/runner/sentieon_runner.py $JSON_FILE) >> $LOG_FILE
echo End: $(date +'%m-%d-%Y-%H-%M-%S') >> $LOG_FILE
