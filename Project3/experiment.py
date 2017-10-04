import sp_approx
import sp_exact_3
import msa_sp_score_3k
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import util


if __name__ == "__main__":
	args = sys.argv
	score, _ = util.read_score_matrix_and_alphabet(args[1])
	gap_cost = int(args[2])
	fileFolder = args[3]
	
	approxScores = []
	exactScores = []
	i = 0
	files = os.listdir(fileFolder)
	files.sort(key=lambda x: int(x.split("_")[1]))
	inputlength = []
	for file in files:
		splitted = file.split("_")

		inputlength.append(int(splitted[1]))
		i += 1
		print("Processing file: {}".format(file))
		sequences = [s.replace(" ", "") for s in util.read_fasta_file(fileFolder + "/" + file).values()]
		
		approxScores.append(sp_approx.compute_score(score, gap_cost, sequences))
		exactScores.append(sp_exact_3.compute_score(score, gap_cost, sequences))
	#print(exactScores)
	#print(approxScores)
		

	with open("results.txt", "w+") as fp:
		for i in range(0,len(exactScores)):
			fp.write(str(inputlength[i]) + "," + str((float(approxScores[i]) / float(exactScores[i]))) + "\n")
			
		
		
