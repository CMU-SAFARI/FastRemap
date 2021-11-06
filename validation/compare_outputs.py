import sys

f1 = open(sys.argv[1]) 
f2 = open(sys.argv[2])

l1 = f1.readline() 
l2 = f2.readline() 

# only care about @SQ for now. 
f1_head = dict() 
f2_head = dict() 

# add unmapped later. 

f1_linenum = 1 
f2_linenum = 1
while True: 
    
    if l1[0] == "@": 
        if "@SQ SN:" in l1: 
            f1_head[l1.split()[1]] = l1.split()[2] 
        l1 = f1.readline() 
        f1_linenum += 1 
    elif l2[0] == "@": 
        if "@SQ SN:" in l2: 
            f2_head[l2.split()[1]] = l2.split()[2] 
        l2 = f2.readline() 
        f2_linenum += 1 
    else: 
        # we are in the body now. 
        # first check the header. 
        for el in f1_head: 
            if el not in f2_head: 
                print("HEADER DIFF: f1 " + str(el) + " -> " + str(f1_head[el])  + "   f2 " + str(el) + " not found") 
            elif f1_head[el] != f2_head[el]: 
                print("HEADER DIFF: f1 " + str(el) + " -> " + str(f1_head[el])  + "   f2 " + str(el) + " -> " + str(f2_head[el]))
        for el in f2_head: 
            if el not in f1_head: 
                print("HEADER DIFF: f1 " + str(el) + " -> " + str(f1_head[el])  + "   f2 " + str(el) + " not found") 
            elif f2_head[el] != f1_head[el]: 
                print("HEADER DIFF: f1 " + str(el) + " -> " + str(f1_head[el])  + "   f2 " + str(el) + " -> " + str(f2_head[el]))
            
        l1_cols_1_11 = l1.split()[0:11] 
        l2_cols_1_11 = l2.split()[0:11] 

        l1_tags = l1.split()[11:]
        l2_tags = l2.split()[11:] 

        exit_code = 0

        for el in l1_tags: 
            if el not in l2_tags: 
                print("TAG DIFF: " + el + " not found in f2") 
                print("f1 line num: " + str(f1_linenum) + "   f2 line num: " + str(f2_linenum))  
                exit_code = 1 

        for el in l2_tags: 
            if el not in l1_tags: 
                print("TAG DIFF: " + el + " not found in f1") 
                print("f1 line num: " + str(f1_linenum) + "   f2 line num: " + str(f2_linenum)) 
                exit_code = 1 
                
        if (exit_code): 
            exit(1) 

        
        for idx,el in enumerate(l1_cols_1_11): 
            if idx == 6: 
                if (l1_cols_1_11[idx] == "=" and (l2_cols_1_11[idx] != l2_cols_1_11[2])) or (l2_cols_1_11[idx] == "=" and (l1_cols_1_11[idx] != l1_cols_1_11[2])): 
                    # RNEXT == "=" seqan3 problem? 
                    print("RNEXT and RNAME DO NOT MATCH")
                    exit_code = 1 
                    break 
                continue 
                    
            if el != l2_cols_1_11[idx]: 
                print("ENTRY DIFF: " + str(el) + " != " + str(l2_cols_1_11[idx])) 
                exit_code = 1 
        
        if exit_code: 
            print("f1 linenum: " + str(f1_linenum) + " --  " + " ".join(l1_cols_1_11)) 
            print("f2 linenum: " + str(f2_linenum) + " --  " + " ".join(l2_cols_1_11)) 
            exit(1)  

        l1 = f1.readline() 
        l2 = f2.readline() 

    if not l1 or not l2: 
        break 

print("Found 0 Errors") 
