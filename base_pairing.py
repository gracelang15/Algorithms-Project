def base_pairing(base_a, base_b):
    if base_a == 'A' or base_a == 'a':
        if base_b == 'U' or base_b == 'u':
            return 1
        else:
            return 0
    elif base_a == 'U' or base_a == 'u':
        if base_b == 'A' or base_b == 'a':
            return 1
        else:
            return 0
    elif base_a == 'G' or base_a == 'g':
        if base_b == 'C' or base_b == 'c':
            return 1
        else:
            return 0
    elif base_a == 'C' or base_a == 'c':
        if base_b == 'G' or base_b == 'g':
            return 1
        else:
            return 0
    #else:
        #return "Wrong input!"
