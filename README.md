# lrf-ms
Tool for computing matching statistics that uses the longest repeated factor (LRF) heuristic.

## Installation

```sh
git clone https://github.com/fmasillo/lrf-ms.git
cd lrf-ms
make
```

## Usage: 

This tool works on fasta files.

```sh
./lrf_ms -r reference_file -c collection_file [options]
  Options: 
    -o <output file> : basename for the output files (default: reference_collection)
    -b <buffer size> : size of the buffer for holding the temporary MS entries in MB (default: 100)
    -t <number of threads> : number of threads to use (default: 8)
    -h : print this help
```


## Citation:

Conference paper:

Zsuzsanna Lipták, Martina Lucà, Francesco Masillo, and Simon J. Puglisi. Fast Matching Statistics for Sets of Long Similar Strings. In Proc. of the Prague Stringology Conference 2024, pages 3--15. Czech Technical University in Prague, Faculty of Information Technology, Department of Theoretical Computer Science, 2024: https://www.stringology.org/event/2024/p02.html

```
@inproceedings{DBLP:conf/stringology/LiptakLMP24,
  author       = {Zsuzsanna Lipt{\'{a}}k and
                  Martina Luca and
                  Francesco Masillo and
                  Simon J. Puglisi},
  editor       = {Jan Holub and
                  Jan Zd{\'{a}}rek},
  title        = {Fast Matching Statistics for Sets of Long Similar Strings},
  booktitle    = {Prague Stringology Conference 2024, Prague, Czech Republic, August
                  26-27, 2024},
  pages        = {3--15},
  publisher    = {Czech Technical University in Prague, Faculty of Information Technology,
                  Department of Theoretical Computer Science},
  year         = {2024},
}
```
