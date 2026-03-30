# Define the File ID from your link
FILE_ID="1FUbW9x23tv-LorTSdOxdXVBn0JRxgRRh"
FILE_NAME="SpCas9_2019_150_sequence.csv"

# 1. Request the download and capture the confirmation code from the cookie
CONFIRM=$(curl -sc /tmp/gcookie "https://drive.google.com/uc?export=download&id=${FILE_ID}" | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1/p')

# 2. Download the file using that confirmation code
curl -Lb /tmp/gcookie "https://drive.google.com/uc?export=download&confirm=${CONFIRM}&id=${FILE_ID}" -o "${FILE_NAME}"

# 3. Clean up the cookie file
rm /tmp/gcookie
