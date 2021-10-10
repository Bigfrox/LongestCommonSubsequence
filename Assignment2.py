'''
Bio Computing Assignment 2, Longest Common Subsequence Finding
2016253072
명수환(Myeong Suhwan)
'''

import re
import random
from datetime import datetime
import string
import sys

def getDataFromFile(filename):
    data1 = ""
    data2 = ""
    try:
        with open(filename, 'r') as file:
            line = None
            line = file.readline() # comment
            if not line: # * input_file is empty
                print("[-] No protein sequence . . .")
                exit(1)
            
            if(line[0] == '>'): # * must accept only first FASTA format.
                print("1st FASTA : [Comment]")
                print(line)
            else:
                print("[-] No correct format . . .")
            while line != '':
                line = file.readline()
                if line:
                    if(line[0] == '>'): # * accept 2nd FASTA format.
                        print("2nd FASTA : [Comment]")
                        print(line)
                        while line != '':
                            
                            line = file.readline()
                            print(line)
                            if line:
                                if(line[0] == '>'): # * ignore 3rd and after FASTA format.
                                    print("[+] 3rd FASTA is ignored . . .")
                                    return data1, data2
                                else:
                                    data2 += line.strip('\n')
                    else:
                        data1 += line.strip('\n')
                else:
                    print("[-] Need one more sequence.")
                    exit(1)
                
    except FileNotFoundError:
        print("[-] No input file . . .")
        
    return data1, data2


def getRandomData(num):
    alphabet = string.ascii_letters
    alphabet = list(alphabet)
    
    protein = alphabet

    data_rand = ""
    
    if num < 1:
        print("[+] Not generated random sequences . . .")
        return data_rand
    if num > 500:
        print("[-] About Maximum recursion depth Error . . ")
        print("[-] This program limit the number of randomly generated data below 500 . . ")
        print("[-] Thus, The number of randomly generated data is set to 500 . . ")
        num = 500
        
    
    
    while num:
        data_rand += protein[random.randint(0,51)]
        num -= 1
    
    return data_rand

def processing(data):
    data = data.replace(" ", "")
    data = data.upper()

    return data

def FindLCSUsingDP(data1, data2):
    
    #! 'i' is index for data 1
    #! 'j' is index for data 2
    m = len(data1)+1
    n = len(data2)+1
    
    score = [[0 for col in range(m)] for row in range(n)]
    
    for i in range(m):
        score[i][0] = 0 # * fill these blocks with 0
    for j in range(n):
        score[0][j] = 0 # * fill these blocks with 0
    for i in range(1,m):
        for j in range(1,n):
            if data1[i-1] == data2[j-1]:
                
                score[i][j] = max(score[i-1][j],score[i][j-1],score[i-1][j-1]+1)
                
            else:
                score[i][j] = max(score[i-1][j],score[i][j-1])
    
    #!DEBUG
    # for i in range(m):
    #     print(score[i])

    subsequence = []
    count = score[m-1][n-1] # * length of LCS
    print("[+] Length of LCS : ",count)
    m -=1
    n -=1
    print("[*] Start the Backtracking . . .")
    BackTracking(data1,score,m,n,subsequence)
    print("[*] End of the Backtracking . . .")
    return subsequence

def BackTracking(data,score, row, col,subsequence):
    
    m = row
    n = col

    if m <= 0 or n <= 0:
        return subsequence

    direction = max(score[m][n-1], score[m-1][n], score[m-1][n-1])
    min_direction = min(score[m][n-1], score[m-1][n], score[m-1][n-1])
    
    if min_direction == score[m][n]: #* No diagonal.
        print("Case - All the same",end = " -> ")
        n -= 1
        print("Now Position : ","(",m ,",",n,")")
        subsequence = BackTracking(data,score,m,n,subsequence)
    elif direction == score[m-1][n-1]: #* Diagonal first
        print("Case - Move diagonally",end = " -> ")
        m -= 1
        n -= 1
        print("Now Position : ", "(",m ,",",n,")")
        subsequence.insert(0,data[m])
        print("[Backtracking] data: ", data[m],"is added.")
        subsequence = BackTracking(data,score,m,n,subsequence)
    elif direction == score[m][n-1]: #* horizontal second
        print("Case - Move horizontally",end = " -> ")
        n -= 1
        print("Now Position : ","(",m ,",",n,")")
        subsequence = BackTracking(data, score,m,n,subsequence)
    else: #* vertical last
        print("Case - Move vertically",end = " -> ")
        m -= 1
        print("Now Position : ","(",m ,",",n,")")
        subsequence = BackTracking(data,score,m,n,subsequence)
    
    if score[m][n] == 0:
        return subsequence



def main():
    if len(sys.argv) != 2:
        print("No input file.")
        print("<Usage> assignment2.py input_filename.txt")
        return -1;
    
    input_filename = sys.argv[1]
    #input_filename = 'test.txt'
    #input_filename = 'nothing.txt'
    #input_filename = '1fasta.txt'
    #input_filename = '3fasta.txt'
    #input_filename = 'no_protein.txt'

    
    data_not = re.compile(r'[^a-zA-Z]')
    data1 = ""
    data2 = ""
    
    data1,data2 = getDataFromFile(input_filename)
    data1 = processing(data1)
    data2 = processing(data2)

    if not data1 or not data2:
        print("No input . . .")
        exit(1)

    no_protein = data_not.findall(data1+data2)
    if no_protein:
        print("No Protein Sequences . . .")
        print(no_protein)
        exit(1)


    num = int(input("Enter the number of Random Data you wanna make : ")) # * Generate Protein Sequences randomly.
    data1 += getRandomData(num)
    data2 += getRandomData(num)

    start_time = datetime.now()
    LCS = FindLCSUsingDP(data1,data2)
    LCS = "".join(LCS)
    #print("Data1 : ",data1)
    #print("Data2 : ",data2)
    print("\n")
    print("[*] Longest Common Subsequence : ",LCS)

    #print(datetime.now())
    print("[+] Time Elapsed : ", datetime.now() - start_time, "microseconds")

if __name__ == '__main__':
    main()