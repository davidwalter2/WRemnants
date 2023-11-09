. ./setup.sh
if [ -e "~/.keytab ${USER}@CERN.CH" ]; then
    echo "Initialize kerberos"
    kinit -kt ~/.keytab ${USER}@CERN.CH
fi
python3 $@
