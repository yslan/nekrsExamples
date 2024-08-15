: ${VERBOSE:=0}
: ${DEBUG:=0}

sess=$1
newfldr=$2
tag=$3

err=0
if [ -z $sess ]; then
  echo "ERROR: <case.sess> not specified"
  err=1
fi
if [ ! -f $sess ]; then
  echo "ERROR: $sess not found"
  err=1
fi
if [ -z $newfldr ]; then
  echo "ERROR: <newfldr> not specified"
  err=2
fi

if [ ! $err -eq 0 ]; then
  echo "exiting ... "$err
  echo -e "\nUsage: [VERBOSE=1] [DEBUG=1] $0 <case.sess> <newfldr> [<tag>]"
  echo -e "\nThis script is for NekNek. It creates new folder, move *0.f0* and copy logfile there\n"

  echo "Variables:"
  echo "    case.sess: .sess file of NekNek case"
  echo "    newfldr: new folder name"
  echo "    tag: tag for logfile"
  echo "    VERBOSE: print what is done"
  echo "    DEBUG: dry-run, only prints command without action"
  exit 1
fi

if [ ! -z $tag ]; then
  tag="_"$tag
fi

if [ $DEBUG -eq 1 ]; then
  VERBOSE=1
fi

while read p; do  # chk neknek subfolders
  # echo $p                   # inlet/inlet:1;
  cname_sess="${p%:*}"        # inlet/inlet
  tmp="${p#*:}"               # 1;
  ntasks_sess="${tmp%;}"      # 1

  if [ -z "$cname_sess" ] || [ -z "$ntasks_sess" ]; then
    echo "ERROR: $case has wrong format at line"
    echo $p
    exit 1
  fi

  fldr="${cname_sess%/*}"     # inlet
  cname="${cname_sess#*/}"    # inlet

  if [ $VERBOSE -eq 1 ]; then
    echo "mkdir -p $newfldr"
    echo "mv $fldr/logfile $newfldr/logfile_$cname$tag"
    echo "mv "$fldr"/"$cname"0.f0* "$newfldr"/"
    echo "mv "$fldr"/"$cname".nek5000 "$newfldr"/"
  fi

  if [ ! $DEBUG -eq 1 ]; then
    mkdir -p $newfldr
    mv $fldr/logfile $newfldr"/logfile_"$cname$tag
    mv `echo $fldr/$cname"0.f0*"` $newfldr"/"
    mv $fldr"/"$cname".nek5000" $newfldr"/"
  fi

done < $sess
