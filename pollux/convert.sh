#!/bin/sh

for file in *.spec
do
  awk '
    NR % 20 == 0 {
      lambda = $1 / 10;
      R = $3;
      if (455 <= lambda && lambda <= 515)
        printf("%.3f %.3f\n", lambda, R);
    }
  ' $file >$file-converted
done
