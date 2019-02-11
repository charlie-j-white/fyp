
# prepare date-time string

string=$(date +"%a")
day=${string^^}
string=$(date +"%d-%b-%Y")
date=${string^^}
string=$(date +"%T")
time=${string^^}

# list of files to add, combine with .gitignore

git add gitsync.sh
git add README.md
git add .gitignore
git add ./frances
git add ./ernie

# rest of the sync commands

git commit -m 'auto sync:'\ $day\ $date\ $time
#git ls-tree -r master --name-only
git push -u origin master





