# Make files read-only
find $root_folder -type f -exec chmod 444 {} \;

# Make folders read/execute-only
find $root_folder -type d -exec chmod 555 {} \;
