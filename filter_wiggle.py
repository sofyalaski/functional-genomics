import sys

def filter_chr(input_file, output_file):

    fileR = open(input_file,'r')
    fileW = open(output_file,'w')

    chromosomes = [str(i) for i in list((range(1,20)))]
    chromosomes.extend(["X","Y"])
    first_line = ["variableStep chrom="+i+'\n' for i in chromosomes]
    first_line_corrected = ["variableStep chrom=chr"+i+'\n' for i in chromosomes]

    for line in fileR:
        if line[0] == 'v':
            if line in first_line:
                #print("match with UCSC lib")
                fileW.write(first_line_corrected[first_line.index(line)])
                for line in fileR:
                    if line[0] !="v":
                        fileW.write(line)
                    else:
                        break
            #else:
                #print(line,"is not in UCSC ")



    fileR.close()
    fileW.close()

def main():
    if (len(sys.argv) != 3):
        print( "two args required")
        exit(-1)
    filter_chr(sys.argv[1],sys.argv[2])


if __name__ == '__main__':
    main()
