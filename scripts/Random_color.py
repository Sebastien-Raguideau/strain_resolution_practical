import random
import sys
import argparse



def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("clusterFile", help="cluster file")
    
    args = parser.parse_args()

    #import ipdb; ipdb.set_trace()    

    bins = set()
    binMap = {}

    bFirst = True
    with open(args.clusterFile,'r') as f:
        for line in f:
            if bFirst:
                bFirst = False
            else:
                line = line.rstrip()

                toks = line.split(',')
                clust = int(toks[1])
                unitig = toks[0]
                bins.add(clust)
            
                binMap[unitig] = clust
            

    number_of_colors = len(bins)

    color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
             for i in range(number_of_colors)]

    binList = list(bins)
    binColor = {}
    
    for b,binL in enumerate(binList):
        binColor[binL] = color[b]
   
    print('NODE,COLOUR') 
    for (unitig,binL) in binMap.items():
        print(unitig + ',' + str(binColor[binL]))

    
if __name__ == "__main__":
    main(sys.argv[1:])
