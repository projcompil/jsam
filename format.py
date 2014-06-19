from sys import argv

print(argv)
[nf, nof, p, q] = argv[1:5]
#p = int(p)
#q = int(q)

with open(nf, "r") as f :
    with open(nof, "a") as outf :
        for (i, line) in enumerate(f) :
            if i > 1 :
                outf.write(line.rstrip() + ",%s,%s\n" % (p,q))
            #elif i == 1 :
            #    outf.write(line.rstrip() + ",fsum,fcorrupt\n")
