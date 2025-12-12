import argparse


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		prog = 'make_ColData ',
		description = "Generate the colData and parse other annotationdb data")
	parser.add_argument('-r',default = [1,2], help = "radius list")
	parser.add_argument('-d',default = [1020,3], help = "dimension list")
	# parser.add_argument('-b', help = "bdb file")
	# parser.add_argument('-j',help = 'jump cp')
	# parser.add_argument('-l', help = "lincs")

	args = parser.parse_args()
	print(args.r)
	print(type(args.r))