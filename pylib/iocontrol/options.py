import sys

def get_options(argvs, chars, datatypes, defaults):
    values=defaults
    for argv in argvs[1::2]:
        if argv in chars:
            datatype=datatypes[chars.index(argv)]
            if datatype=='int':
                values[chars.index(argv)]=int(argvs[argvs.index(argv)+1])
            elif datatype=='float':
                values[chars.index(argv)]=float(argvs[argvs.index(argv)+1])
            elif datatype=='str':
                values[chars.index(argv)]=argvs[argvs.index(argv)+1]
            else:
                print('A data type given in the datatype array, '+datatype+' is not valid.')
                sys.exit()
        else:
            print('The option given in the command line, '+argv+' is not valid.')
            options='['
            for op in chars:
                options+=op
                if op!=chars[-1]:
                    options+=', '
            options+=']'
            print('The possible options are '+options+'.')
            sys.exit()

    return values



