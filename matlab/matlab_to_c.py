# this should work for one char variables
eq = input("Enter string\n")
operators = ['+', '-', '/', ' ', ')', '*', '^', '(']
add_i = len("pow(,x)") -1
out = []
pow_idx = [];
i = 0
while i < len(eq):
    c = eq[ i ]
    # print(c)
    if c == '^':
        load = True
        left_b = 0
        idx_exp = 0 
        idx_minus = 0
        buf = []
        while load:
            idx_exp += 1
            buf.append(eq[ i + idx_exp ])
            if eq[i + idx_exp] == '(':
                left_b += 1
            elif eq[i + idx_exp] == ')':
                left_b -= 1
            if left_b == 0:
                load = False
        if buf[0] == '(' and buf[idx_exp -1] == ')':
            buf = buf[1: idx_exp -1 ]
            idx_minus += 2
        

        power = ''.join(buf)
        load = True
        right_b = 0
        idx = 0 
        buf = []
        while load:
            idx +=1
            buf.insert(0, eq[i - idx])
            if eq[i - idx] == ')':
                right_b += 1
            elif eq[i - idx] == '(':
                right_b -= 1
            if right_b == 0:
                load = False
                if len(buf) == 1:
                    while eq[i - idx - 1] not in operators:
                        buf.insert(0, eq[i - idx -1])
                        idx +=1

        if buf[0] == '(' and buf[idx-1] == ')':
            buf = buf[1: idx -1]
            idx_minus += 2
        buf = ''.join(buf)
        # eq = eq[0: i-idx] + 'pow(' +   ')' + eq[i+2]
        if power == '1/2':
            eq = eq[0: i-idx] + 'sqrt(' + buf + ')' + eq[i + idx_exp + 1:len(eq)]
            i += 4
        else:
            eq = eq[0: i-idx] + 'pow(' + buf + ',' + power + ')' + eq[i + idx_exp + 1:len(eq)]
            i +=7
        i -= idx_minus
    i += 1 

print(eq)

                   

