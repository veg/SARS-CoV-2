#!/usr/bin/zsh

source $HOME/.zshrc

echo "Adding newlines where headers are"
sed -i 's/\(\w\)>/\1\n>/' $1
echo "Converting DOS newlines to Unix"
dos2unix $1
echo "Uploading metadata to MongoDB"
node submit-to-mongo.js $2
echo "Uploading sequences to MongoDB"
python3 /data/shares/new-gisaid/update_with_sequence.py -i $1
