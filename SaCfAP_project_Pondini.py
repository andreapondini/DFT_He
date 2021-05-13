R_MIN, R_MAX, SAMPLES, IOUNIT = 0, 50 , 4049, 10

#only closed shell s elements are supported (He & Be)
atom_type = " "
while atom_type != "He" and atom_type != "Be":
    atom_type = input("Choose either the He or Be atom: ")
    if atom_type != "He" and atom_type != "Be":
        print("You must type He or Be")
if atom_type == "He":
    #He
    NUCLEAR_CHARGE = 2
    N_ELECTRONS = 2
elif atom_type == "Be":
    #Be
    NUCLEAR_CHARGE = 4
    N_ELECTRONS = 4
    
