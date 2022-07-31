# This is the github tutorial
# created 2021 Dec, modified 2022 Feb

# github keychain
# GitHub keychain 
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


### the whole pipeline
cd into current directory you wanna upload

	git init

add a git URL as an alias

	git remote add origin https://github.com/XiaoqiongGu/Gu_2021_Augmentin16S.git

## make changes
browse and inspect the evolution of project files

	git log -1 -p
	git reflog

show modified files in working directory,staged for your next commit

	git status

	git diff

add or update a file as it looks now to your next commit

	git add .
	git add -A  # Flags and parameters in git are case sensitive, and the -A flag can stage all my unstaged changes.

commit your staged content as a new committ snapshot

	git commit -m 'descriptive message'


push the committed change back to github

	git push
	git push --set-upstream origin main # the first time you push it back to your source control, you need to include the

clone the repo

	git clone 

if you already clone the repo and wanna to retrive new changes

	git pull --rebase

checkout a new branch (-b to create it). You can also use checkout to move between branches on your local system.

	git checkout -b new-branch

or if you wanted to move from a branch you created back to the master branch

	git checkout master

nano .gitignore

	# OSX files
	.DS_Store
	
	# VS Code files
	.vscode
	
	# Jupyter Notebook files
	.ipynb_checkpoints


	git stash 

	git stash temporarily shelves (or stashes) changes you've made to your working copy so you can work on something else, and then come back and re-apply them later on. Stashing is handy if you need to quickly switch context and work on something else, but you're mid-way through a code change and aren't quite ready to commit.

	gitk 

	git reset

	git reset -HEAD^ 
	git reset commitid -hard  ### pick 

	git cherry.pick

	git clean -f # very fatal, and can be very dangerous

	git rm vs rm 

	what is the difference between git pull and git clone
		so depending on the location 
		git pull is widely used often

where can you get more information?

# git revert back

[how do i revert a git repo to a previous commit](https://stackoverflow.com/questions/4114095/how-do-i-revert-a-git-repository-to-a-previous-commit)


# Git 20190416
任何类型文件进行版本控制  
1. 本地版本控制系统  
2. 集中化的版本控制系统  
3. 分布式版本控制系统  

## Git usage

# 20190523
1. make a new repo on the github website, can choose public and private repo
2. follow instructions on the website,
  echo "# underworldvirome" >> README.md
  git init
  git add README.md
  git commit -m "first commit"
  git remote add origin https://github.com/XiaoqiongGu/underworldvirome.git
  git push -u origin master
3. when files added in the folder
  git add.
  git commit
  git push origin master
  git status
