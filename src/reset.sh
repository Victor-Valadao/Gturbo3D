#!/bin/bash

rm ./out*
rm ./slurm*
rm ./fields/*
mv ./files/seed.000 ./seed.000
rm ./files/*
mv ./seed.000 ./files/seed.000
rm ./curframe.dat
cp ../src/startframe.dat curframe.dat
