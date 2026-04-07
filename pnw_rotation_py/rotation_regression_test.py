import csv
import sys
import numpy as np
from geo_helper import GeoHelper as gh

def main():
    gh.clamp(1,2,3)
    source = input("Enter source: 1 - Wells simpson, 2 - large set, 3 - quadratic test")
    if source == "1":
        print("parameters from csv from Wells Simpson paper (edited)")
        file_path = "data/testData.csv"
        print(f"Selected file: {file_path}")
        with open(file_path, mode='r', newline='', encoding='utf-8') as file:
            reader = csv.DictReader(file)
            #print(reader.fieldnames)
            TR = []
            for row in reader:
                if row["R"] != "" and row["T"] != "":
                    TR.append([float(row["R"]),float(row["T"])])

    if source == "2":
        print("test parameters from csv using published paleo")
        file_path = "data/Unit,Sub-unitSite,Age_Ma,Rotation_d.csv"
        print(f"Selected file: {file_path}")
        with open(file_path, mode='r', newline='', encoding='utf-8') as file:
            reader = csv.DictReader(file)
            #print(reader.fieldnames)
            TR = []
            for row in reader:
                if row["Rotation_deg"] != "" and row["Age_Ma"] != "":
                    TR.append([float(row["Rotation_deg"]), float(row["Age_Ma"])])

    if source == "3":
        print("test parameters from quadratic fit to graph:")
        a = -0.0088 # from quadratic fit
        b = 0.4766  # from quadratic fit
        w_0 = 1.19  # from Wells - Simpson Fig3B
        print("a = ", a)
        print("b = ", b)
        print("w_o = ", w_0)
        TR = createTestDataSet(a, b, w_0, 50, 10)

    x = calculateRateCoeff(TR, normalize=True)
    print("x: ", x)

def calculateRateCoeff(trArray, normalize = True):
    T1 = T2 = T3 = T4 = T5 = T6 = RT1 = RT2 = RT3 = 0.0
    
    for TR in trArray:
        T = TR[1]
        R = TR[0]
        print (TR)
        T1 += T
        T2 += pow(T, 2)
        T3 += pow(T, 3)
        T4 += pow(T, 4)
        T5 += pow(T, 5)
        T6 += pow(T, 6)

        RT1 += R * T
        RT2 += R * T * T
        RT3 += R * T * T * T

    A = np.array([
        [T6 / 3, T5 / 2, T4],
        [T5 / 3, T4 / 2, T3],
        [T4 / 3, T3 / 2, T2]
    ])
    #print("A: ", A)

    b = np.array([RT3, RT2, RT1])
    #print("b: ", b)

    x = np.linalg.solve(A, b)
    if (normalize):
        x = x / x[2]
    return x

def createTestDataSet(a, b, w_0, tMax, N):
    RT = []
    deltaT = tMax / N
    for i in range(N):
        T = i * deltaT
        R = (1.0/3.0 * a * pow(T,3) + 0.5 * b * pow(T,2) + T) * w_0
        RT.append([R, T])
    return RT

if __name__ == "__main__":
    # The return value of main() is passed to sys.exit()
    # to provide a proper exit code to the operating system.
    sys.exit(main())