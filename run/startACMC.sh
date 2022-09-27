#!/bin/bash
cp Backup_InitialData InitialData
sed -i  '27 c\1 Create_ID' Config
../src/ACMC_Data
cp Created_InitialData InitialData
sed -i  '27 c\0	Create_ID' Config
../src/ACMC_Data
