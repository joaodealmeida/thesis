#!/usr/bin/python

import csv
from sys import argv

def main():
	print '>>> Opening file...'

	with open(argv[1]) as tsv:
		print tsv
		print

		print '>>> Parsing file...'
		for line in csv.reader(tsv, delimiter='\t'):
			print line
		print

	print '>>> Everything done!'


if __name__ == '__main__': main()
