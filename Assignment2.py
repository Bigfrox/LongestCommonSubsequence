'''
Bio Computing Assignment 2, Longest Common Subsequence Finding
2016253072
명수환(Myeong Suhwan)

'''

import re
import random
from datetime import datetime
import string

def getDataFromFile(filename):
    data1 = ""
    data2 = ""
    try:
        with open(filename, 'r') as file:
            line = None
            line = file.readline() # comment
            if not line:
                print("No input data . . .")
                exit(1)
            
            if(line[0] == '>'): # * must accept only first FASTA format.
                print("1st FASTA : [Comment]")
                print(line)
            else:
                print("No correct format . . .")
            while line != '':
                line = file.readline()
                if line:
                    if(line[0] == '>'): # * accept 2nd FASTA format.
                        print("2nd FASTA : [Comment]")
                        print(line)
                        while line != '':
                            line = file.readline()
                            if line:
                                if(line[0] == '>'): # * ignore 3rd and after FASTA format.
                                    break
                                else:
                                    data2 += line.strip('\n')
                    else:
                        data1 += line.strip('\n')
                
    except FileNotFoundError:
        print("No input file . . .")
        
    return data1, data2


def getRandomData(num):
    alphabet = string.ascii_letters
    alphabet = list(alphabet)
    
    protein = alphabet
    if num < 1:
        print("Not generated random sequences . . .")
        
    data_rand = ""
    
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
        score[i][0] = 0
    for j in range(n):
        score[0][j] = 0
    for i in range(1,m):
        for j in range(1,n):
            if data1[i-1] == data2[j-1]:
                
                score[i][j] = max(score[i-1][j],score[i][j-1],score[i-1][j-1]+1)
                
            else:
                score[i][j] = max(score[i-1][j],score[i][j-1])
    
    #!DEBUG
    for i in range(m):
        print(score[i])

    subsequence = []
    count = score[m-1][n-1] # * length of LCS
    print("count : ",count)
    m -=1
    n -=1
    BackTracking(data1,score,m,n,subsequence)
    return subsequence

def BackTracking(data,score, row, col,subsequence):
    m = row
    n = col

    if m <= 0 or n <= 0:
        return subsequence

    direction = max(score[m][n-1], score[m-1][n], score[m-1][n-1])
    min_direction = min(score[m][n-1], score[m-1][n], score[m-1][n-1])
    
    if min_direction == score[m][n]: #* No diagonal.
        print("다같은경우")
        n -= 1
        print("현재 위치 : ", m ,n)
        subsequence = BackTracking(data,score,m,n,subsequence)
    elif direction == score[m-1][n-1]: #* Diagonal first
        print("대각선 이동")
        m -= 1
        n -= 1
        print("현재 위치 : ", m ,n)
        subsequence.insert(0,data[m])
        print("data: ", data[m],"을 추가하였습니다.")
        subsequence = BackTracking(data,score,m,n,subsequence)
    elif direction == score[m][n-1]: #* horizontal second
        print("수평 이동")
        print("현재 위치 : ", m ,n)
        n -= 1
        subsequence = BackTracking(data, score,m,n,subsequence)
    else: #* vertical last
        print("수직 이동")
        print("현재 위치 : ", m ,n)
        m -= 1
        subsequence = BackTracking(data,score,m,n,subsequence)
    
    if score[m][n] == 0:
        return subsequence



def main():
    

    
    data_not = re.compile(r'[^a-zA-Z]')
    data1 = ""
    data2 = ""
    
    filename = 'test.txt'
    #! filname = argv
    
    data1,data2 = getDataFromFile(filename)
    num = 0 # * Generate num DNA Sequences.
    data1 += getRandomData(num)
    data2 += getRandomData(num)
    
    
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

    #! for debug
    #data1 = "ATGTTAT"
    #data2 = "ATCGTAC"

    print("[Data1]\n",data1)
    print("\n")
    print("[Data2]\n",data2)
    #print(len(data))

    start_time = datetime.now()
    LCS = FindLCSUsingDP(data1,data2)
    LCS = "".join(LCS)
    print("[Data1]\n",data1)

    print("[Data2]\n",data2)
    print("Longest Common Subsequence : ",LCS)

    #print(datetime.now())
    print("Time Elapsed : ", datetime.now() - start_time, "microseconds")

if __name__ == '__main__':
    main()