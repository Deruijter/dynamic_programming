# Author: Markus de Ruijter

# This class contains methods to perform different kinds of alignment
# Properties should be set when instantiating the class
# The public method align() is then available to perform the alignments

from entities.exchange_matrices import ExchangeMatrices

class AlignmentHandler:
	
	# Private properties
	__coordinates_aligned_localy = []
	__alignment_score = ""
	__alignment_seq1 = ""
	__alignment_seq2 = ""
	__alignment_similarity = ""
	__all_local_alignments_found = False
	__is_first_local_alignment = True

	__traceback_start_row = 0
	__traceback_start_col = 0

	# Initializer for this class. Class properties should be set here.
	def __init__(self, alignment_type, exchange_matrix_name, sequences, penalty, verbose, local_alignment_limiter_type, local_alignment_limiter_value):
		self.__sequence_1 = sequences[0].Sequence
		self.__sequence_2 = sequences[1].Sequence	
		self.__penalty = penalty
		self.__verbose = verbose

		self.__set_alignment_type(alignment_type)
		self.__set_exchange_matrix(exchange_matrix_name)
		self.__set_local_alignment_limiter(local_alignment_limiter_type,local_alignment_limiter_value)

		self.__clear_class_variables()

	# General alignment method, call other methods based on the alignment type
	def align(self):
		# For local alignment we can have multiple alignments
		if self.__local_alignment:
			while not self.__all_local_alignments_found:
				self.__clear_class_variables()

				self.__create_scoring_matrix()
				self.__set_traceback_start_position()
				self.__do_traceback()
				self.__print_results()

				self.__is_first_local_alignment = False
				self.__local_alignment_limit_iterations -= 1
				if self.__local_alignment_limit_iterations <= 0:
					break
		else:
			self.__create_scoring_matrix()
			self.__set_traceback_start_position()
			self.__do_traceback()
			self.__print_results()

	# Clear class variables so we can re-run the process if needed
	def __clear_class_variables(self):
		self.__alignment_score = ""
		self.__alignment_seq1 = ""
		self.__alignment_seq2 = ""
		self.__alignment_similarity = ""
		self.__all_local_alignments_found = False
		self.__traceback_start_row = 0
		self.__traceback_start_col = 0

	# Print the scoring matrix, sequence alignment and the alignment score
	def __print_results(self):
		if self.__local_alignment:			
			if self.__alignment_seq1 == "" and self.__is_first_local_alignment:
				print("No proper local alignment could be made")
				return
			elif len(self.__alignment_seq1) < self.__local_alignment_limit_length:
				return
			elif self.__all_local_alignments_found:
				return

		self.__print_score_matrix()
		print(self.__alignment_seq1)
		print(self.__alignment_similarity)
		print(self.__alignment_seq2)
		print("Alignment score: {0}".format(self.__score_matrix[self.__traceback_start_row][self.__traceback_start_col]))
		print("")

	# Do alignment based on alignmenttype and return a result matrix
	def __create_scoring_matrix(self):
		range_seq1 = range(0, len(self.__sequence_1)+1)	    
		range_seq2 = range(0, len(self.__sequence_2)+1)
		self.__score_matrix = [[0 for col in range_seq2] for row in range_seq1]

		self.__add_matrix_starting_gaps()
		self.__fill_score_matrix()

	# Set the position to start the traceback from, depending on the alignment type
	def __set_traceback_start_position(self):
		rows = range(len(self.__score_matrix))
		columns = range(len(self.__score_matrix[0]))
		last_row = len(self.__score_matrix)-1
		last_col = len(self.__score_matrix[0])-1
		highest_value = 0;

		# Global alignment just takes the lower right corner of the matrix
		if self.__global_alignment:
			self.__traceback_start_row = last_row
			self.__traceback_start_col = last_col

		# Semi-global takes the highest value on the last row or column
		elif self.__semi_global_alignment:
			for row in rows:
				if self.__score_matrix[row][last_col] >= highest_value:
					highest_value = self.__score_matrix[row][last_col]
					self.__traceback_start_row = row
					self.__traceback_start_col = last_col
			for col in columns:
				if self.__score_matrix[last_row][col] >= highest_value:
					highest_value = self.__score_matrix[last_row][col]
					self.__traceback_start_row = last_row
					self.__traceback_start_col = col

		# Local alignment uses the highest value of the whole matrix
		elif self.__local_alignment:
			for row in rows:
				for col in columns:
					if self.__score_matrix[row][col] >= highest_value:
						highest_value = self.__score_matrix[row][col]
						self.__traceback_start_row = row
						self.__traceback_start_col = col

			# In case we have found all (required) alignments
			if highest_value < self.__local_alignment_limit_score or highest_value == 0:
				self.__all_local_alignments_found = True

	# Perform the part of the traceback that actually checks matches
	def __do_traceback(self):
		if self.__local_alignment and self.__all_local_alignments_found:
			return

		last_row = len(self.__score_matrix)-1
		last_col = len(self.__score_matrix[0])-1

		current_row = self.__traceback_start_row
		current_col = self.__traceback_start_col

		# Handle end gaps for semi-global alignment
		if self.__semi_global_alignment:
			if current_row == last_row:
				amount_of_end_gaps = last_col - current_col
				for gap in range(amount_of_end_gaps):
					self.__alignment_seq1 += "-"
					self.__alignment_seq2 += ("-"+self.__sequence_2)[last_col - gap]
					self.__alignment_similarity += " "
			elif current_col == last_col:
				amount_of_end_gaps = last_row - current_row
				for gap in range(amount_of_end_gaps):
					self.__alignment_seq1 += ("-"+self.__sequence_1)[last_row - gap]
					self.__alignment_seq2 += "-"
					self.__alignment_similarity += " "

		# Handle the actual values
		while current_row > 0 or current_col > 0:
			current_score = self.__score_matrix[current_row][current_col]

			# Because the starting gaps are free for semi-global alignment
			penalty = self.__penalty
			if self.__semi_global_alignment and (current_col == 0 or current_row == 0):
				penalty = 0

			# all the scores of the 3 possible paths (diagonal, horizontal, vertical)
			score_d = self.__score_matrix[current_row-1][current_col-1]
			score_h = self.__score_matrix[current_row][current_col-1]
			score_v = self.__score_matrix[current_row-1][current_col]

			# required scores for the 3 paths in order to be taken
			match_score = self.__exchange_matrix[ord(self.__sequence_1[current_row-1]) - ord("A")][ord(self.__sequence_2[current_col-1]) - ord("A")]
			req_score_d = current_score - match_score
			req_score_h = current_score + penalty
			req_score_v = current_score + penalty

			# We're taking the highroad by default
			if score_v == req_score_v and current_row >0:
				# Store coordinates for future iterations (local alignment)
				self.__coordinates_aligned_localy.append((current_row, current_col))

				self.__alignment_seq1 += self.__sequence_1[current_row-1]
				self.__alignment_seq2 += "-"
				self.__alignment_similarity += " "

				current_row -= 1

				if self.__local_alignment and score_v == 0:
					break

			elif score_d == req_score_d:
				self.__coordinates_aligned_localy.append((current_row, current_col))

				self.__alignment_seq1 += self.__sequence_1[current_row-1]
				self.__alignment_seq2 += self.__sequence_2[current_col-1]
				if self.__sequence_1[current_row-1] == self.__sequence_2[current_col-1]:
					self.__alignment_similarity += "|"
				else:
					self.__alignment_similarity += " "

				current_row -= 1
				current_col -= 1

				if self.__local_alignment and score_d == 0:
					break

			elif score_h == req_score_h and current_col >0:
				self.__coordinates_aligned_localy.append((current_row, current_col))

				self.__alignment_seq1 += "-"
				self.__alignment_seq2 += self.__sequence_2[current_col-1]
				self.__alignment_similarity += " "

				current_col -= 1

				if self.__local_alignment and score_h == 0:
					break

			else:
				print("Error: Traceback cannot be completed properly")
				print("Stuck at {0}, {1}".format(current_row, current_col))
				print("Score v-d-h: {0}, {1}, {2}".format(score_v, score_d, score_h))
				print("Req score v-d-h: {0}, {1}, {2}".format(req_score_v, req_score_d, req_score_h))
				self.__print_score_matrix()
				exit(-1)

		# reverse the strings so we have the originals again
		self.__alignment_seq1 = self.__alignment_seq1[::-1]
		self.__alignment_seq2 = self.__alignment_seq2[::-1]
		self.__alignment_similarity = self.__alignment_similarity[::-1]


	# Print a matrix
	def __print_score_matrix(self):
		if self.__verbose:
			matrix_with_headers = self.__add_matrix_headers()
			print("\n".join([''.join(['{:>4}'.format(item) for item in row]) for row in matrix_with_headers]))

	# Fill the score matrix based on the alignment type chosen
	# This method does not handle the initial gaps
	def __fill_score_matrix(self):
		# Start at 1, the starting gaps are handled by a different method
		rows = range(1, len(self.__score_matrix))
		columns = range(1, len(self.__score_matrix[0]))

		for row in rows:
			for col in columns:
				# Calculate 3 scores (horizontal, vertical, diagonal)
				match_score = self.__exchange_matrix[ord(self.__sequence_1[row-1]) - ord("A")][ord(self.__sequence_2[col-1]) - ord("A")]
				score_d = self.__score_matrix[row-1][col-1] + match_score
				score_h = self.__score_matrix[row][col-1] - self.__penalty
				score_v = self.__score_matrix[row-1][col] - self.__penalty

				if self.__global_alignment or self.__semi_global_alignment:
					self.__score_matrix[row][col] = max([score_h, score_v, score_d])

				elif self.__local_alignment:
					# If we have already aligned this pair, set it to 0 so we won't align it again
					if (row, col) in self.__coordinates_aligned_localy:
						self.__score_matrix[row][col] = 0
						continue

					self.__score_matrix[row][col] = max([score_h, score_v, score_d, 0])


	# Add the headers (sequences) to the matrix.
	# Note: Header must be added after calculations have run
	# Return type: matrix
	def __add_matrix_headers(self):		
		total_rows = range(len(self.__score_matrix) +1)
		total_cols = range(len(self.__score_matrix[0]) +1)
		matrix = [[0 for col in total_cols] for row in total_rows]

		for row in total_rows:
			matrix[row][0] = (" -"+self.__sequence_1)[row]

		for col in total_cols:
			matrix[0][col] = (" -"+self.__sequence_2)[col]

		for row in total_rows:
			for col in total_cols:
				if row > 0 and col > 0:
					matrix[row][col] = self.__score_matrix[row-1][col-1]

		return matrix

	# Add starting gaps to the matrix, depending on the alignment type
	def __add_matrix_starting_gaps(self):
		total_rows = range(1, len(self.__score_matrix))
		total_cols = range(1, len(self.__score_matrix[0]))
 		
		self.__score_matrix[0][0] = 0
		
		for row in total_rows:
			if self.__global_alignment:
				value = self.__score_matrix[row-1][0] - self.__penalty
			elif self.__semi_global_alignment or self.__local_alignment:
				value = 0
			self.__score_matrix[row][0] = value

		for col in total_cols:
			if self.__global_alignment:
				value = self.__score_matrix[0][col-1] - self.__penalty
			elif self.__semi_global_alignment or self.__local_alignment:
				value = 0
			self.__score_matrix[0][col] = value

	# Set limit values for local alignment so we can reduce the amount of results
	def __set_local_alignment_limiter(self, ltype, lvalue):
		if not self.__local_alignment:
			return

		self.__local_alignment_limit_score = 0
		self.__local_alignment_limit_length = 0
		self.__local_alignment_limit_iterations = 999 # Just in case

		if ltype == "score":
			print("Only printing alignments with a score >= {0}".format(lvalue))
			self.__local_alignment_limit_score = lvalue
		elif ltype == "length":
			print("Only printing alignments with a minimum length of {0} chars".format(lvalue))
			self.__local_alignment_limit_length = lvalue
		elif ltype == "iterations":
			print("Printing the top {0} scores".format(lvalue))
			self.__local_alignment_limit_iterations = lvalue

		print("")

	# Set alignment properties so we can use them more efficiently in code
	def __set_alignment_type(self, alignment_type):
		self.__global_alignment = False
		self.__semi_global_alignment = False
		self.__local_alignment = False

		if alignment_type == "global":
			self.__global_alignment = True
		elif alignment_type == "semi-global":
			self.__semi_global_alignment = True
		elif alignment_type == "local":
			self.__local_alignment = True
		else:
			print("this alignment type has not yet been implemented")
			exit(-1)
		return

	# Get the actual exchange matrix we will use
	def __set_exchange_matrix(self, exchange_matrix_name):
		try:
			self.__exchange_matrix = getattr(ExchangeMatrices, exchange_matrix_name)
		except(Exception, e):
			print("unknown exchange matrix", exchange_matrix_name)
			exit(-1)
