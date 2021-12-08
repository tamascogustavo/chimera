# -*- coding: utf-8 -*-
# @Author: gustavotamascohotmail.com
# @Date:   2021-11-30 15:01:26
# @Last Modified by:   tamascogustavo
# @Last Modified time: 2021-12-08 09:09:09

import os 
import subprocess
from sys import argv
import argparse
import sys

def menu():
	'''
	This function creates the script menu
	'''
	parser = argparse.ArgumentParser(description='Add new media to CarveMe medium database.')
	parser.add_argument("--add_media", help="gets a tsv file containing new media description and \
		add to media_db.tsv of CarveMe", default="<new_media.tsv>.")
	args = parser.parse_args()
	return args

def find_file(query):
	'''
	This function find the path to media_db file

	:param query: is the name of the file to be found
	:return all_file_path_matches: a list of all possible files
	'''
	all_file_path_matches = []
	for root, dirs, files in os.walk(r"/"):
		for name in files:
			if name == query:
				hit = os.path.abspath(os.path.join(root, name))
				all_file_path_matches.append(hit)
	return all_file_path_matches

def parse_file(f):
	'''
	This function creates a databese of media

	:param f : is a tsv file
	:return new_media: a list of the components present in the media file provided
	
	'''
	new_media = []
	for index, line in enumerate(f):
		if index>0:
			medium = line.strip().split("\t")[0]
			description = line.strip().split("\t")[1]
			compound = line.strip().split("\t")[2]
			name = line.strip().split("\t")[3]
			item = [medium, description, compound, name]
			new_media.append(item)
	return new_media

def compile_info(carve_db, new_media):
	'''
	This function combine carveme db with your media

	:param carve_db: list of carveme_media_db 
	:param new_media: list of your_media_db

	:return database: an updated carveme media db
	'''
	database = parse_file(carve_db)
	database.extend(new_media)
	return database

def save_database(db_list, path):
	'''
	This function saves the updated media to carveme media db
	
	:param path: path to media_db.tsv from carveme
	:param db_list: is the updated list of media
	'''
	sys.stdout = open(path, "w")
	print("medium\tdescription\tcompound\tname")
	for item in db_list:
		print(f"{item[0]}\t{item[1]}\t{item[2]}\t{item[3]}")

	sys.stdout.close()

def main():
	dir_path = os.getcwd()
	arguments = menu()
	pattern = "/anaconda3/envs/chimera/lib/python3.7/site-packages/carveme/data/input"
	target = "media_db.tsv"
	new_media_file = arguments.add_media
	possible_hits = find_file(target)
	media_file_path = [x for x in possible_hits if pattern in x]
	if len(media_file_path) > 1:
		media_file_path = media_file_path[0]
	else:
		media_file_path = media_file_path

	with open(new_media_file) as n_m_file:
		new_media_info = parse_file(n_m_file)

	with open(media_file_path) as carve_db:
		new_databe = compile_info(carve_db, new_media_info)
	save_database(new_databe, media_file_path)
	



if __name__ == '__main__':
	main()