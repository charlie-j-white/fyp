
# prepare date-time string

time=$(date +"%a %d-%b%Y %T")
string=${time^^}

# list of files to add, combine with .gitignore

git add gitsync.sh
git add README.md
git add ./frances
git add ./ernie

# rest of the sync commands

git commit -m '$string'
#git ls-tree -r master --name-only
git push -u origin master





