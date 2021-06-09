import os

def readSequenceFile(filepath):
    sequence_map = {}
    try:
        open(filepath)
    except NameError:
        print("No input file")
    else:
        if os.path.exists(filepath):
            with open(filepath, 'rb') as f:
                try:
                    content = f.read()
                    content = content.splitlines()
                    count = 0
                    for i in range(len(content)):
                        if (content[i]) and (content[i] is not None) and (len(content[i]) >= 1000):
                            m = str(content[i], 'utf-8')
                            b = True
                            string= 'ATGC'
                            for j in range(len(m)):
                                if m[j] not in string:
                                    # print(m[j])
                                    b = False
                                    break
                            if b == True:
                                sequence_map[count] = m
                                count+=1
                        else:
                            print("Invalid sequence: ", str(i))
                except OSError:
                    print("Error in reading input sequence file.")
    return sequence_map

