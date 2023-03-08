#!/bin/bash

# crontab -e
# 0 2 * * * ~/clean_tmp.sh
# /var/spool/cron/crontabs/

rm /tmp/*.log
