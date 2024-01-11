# unfolding bias test
# 0) Generate card with nominal model and pseudodata sets with different masses 
# 1) Unfold pseudodata with mW set to values of {0, 10, 20, ...}
# 2) Generate card with nominal theory model
# 3) Fit theory model to unfolded pseudo data 
# 4) Compare observed mass pulls with pseudodata value

# hyperparameters
# a) uncertainty on mass when unfolding
# a1) None
# a2) 100MeV
# b) Histmaker type
# b1) mz_wlike
# b2) mw
# c) transverse mass cut (only for mz_wlike possible on real data)
# c1) with cuts
# c2) without cut 

# input arguments: 

HISTMAKER_DIR=$1
COMBINE_DIR=$2

shift
shift
OPTSTAT=""  # other options can be: --massUnc X; --addTauToSignal; --doStatOnly
OPTMASS=""
OPTTAU=""
STATONLY=""
MASS=""
for arg in "$@"; do
  # Check if the current argument is the specific string you're looking for
  if [ "$arg" = "--doStatOnly" ]; then
    echo "Run with stat only."
    STATONLY=_statOnly
    OPTSTAT="--doStatOnly"
  elif [ "$arg" = "--addTauToSignal" ]; then
    echo "Run with stat only."
    OPTTAU="--addTauToSignal"
  elif [ "$arg" = "--massVariation" ]; then
    echo "Run with mass uncertainty."
    MASS="_massVariation"
  elif [ "$MASS" = "_massVariation" ]; then
    echo "Mass uncertainty is ${arg} MeV."
    MASS="${MASS}${arg}"
    OPTMASS="--massVariation ${arg}"
  fi

done

ANALYSIS=WMass # ZMassWLike
LABEL=${ANALYSIS:0:1}

# settings
GEN_VARS="qGen-ptGen" #-absEtaGen"
GEN_AXES=$(echo "$GEN_VARS" | sed 's/-/ /g')
GEN_VARS_STR=$(echo "$GEN_VARS" | sed 's/-/_/g')
echo "GEN_VARS=${GEN_VARS}"
echo "GEN_AXES=${GEN_AXES}"
echo "GEN_VARS_STR=${GEN_VARS_STR}"

FIT_VARS="pt-charge"
FIT_VARS_STR=$(echo "$FIT_VARS" | sed 's/-/_/g')
FIT_VARS_STR="${FIT_VARS_STR}_passIso_passMT"
FAKERATE_VARS="pt charge"
echo "FIT_VARS=${FIT_VARS}"
echo "FIT_VARS_STR=${FIT_VARS_STR}"
echo "FAKERATE_VARS=${FAKERATE_VARS}"

CMSSW_BASE=/home/${USER:0:1}/${USER}/CMSSW_10_6_30/src/


# Do the same thing for multiple histmakers found in the directory
for HISTMAKER_FILE in "$HISTMAKER_DIR"/*; do
    # Check if the item is a file
    if [ ! -f "$HISTMAKER_FILE" ]; then
        echo "Skipping non-file: $HISTMAKER_FILE"
        continue
    fi
    echo "Now at histmaker $HISTMAKER_FILE"

    # suffix from the histmaker
    SUFFIX=$(echo "$HISTMAKER_FILE" | awk -F"Corr_" '{print $2}')
    SUFFIX=${SUFFIX%.hdf5}

    if [ "$SUFFIX" = "" ]; then
        SUFFIX="nominal"
    fi

    echo "suffix ${SUFFIX}"

    COMBINE_OUTDIR=${COMBINE_DIR}/${SUFFIX}

    # 0) Generate cards with different settings with nominal model and pseudodata sets with different masses 
    for massVariation in 0 10 100; do

        OPTMASS="--massVariation $massVariation"

        COMBINE_ANALYSIS_OUTDIR=${COMBINE_OUTDIR}/${ANALYSIS}_${FIT_VARS_STR}_massVar${massVariation}${STATONLY}/
        COMBINE_ANALYSIS_PATH=${COMBINE_ANALYSIS_OUTDIR}/${ANALYSIS}.hdf5

        # 0) Generate cards with nominal model and pseudodata sets with different masses 
        if [ -e $COMBINE_ANALYSIS_PATH ]; then
            echo "The file $COMBINE_ANALYSIS_PATH exists, continue using it."
        else
            echo "The file $COMBINE_ANALYSIS_PATH does not exists, produce it."
            COMMAND="./scripts/ci/run_with_singularity.sh scripts/ci/setup_and_run_python.sh scripts/combine/setupCombine.py \
                -i $HISTMAKER_FILE -o $COMBINE_OUTDIR --hdf5 --sparse --unfolding --ABCD $OPTTAU \
                --genAxes $GEN_AXES --fitvar $FIT_VARS --fakerateAxes $FAKERATE_VARS \
                --pseudoData massWeight$LABEL --pseudoDataAxes massShift --pseudoDataIdxs -1 $OPTSTAT $OPTMASS"
            echo $COMMAND
            eval $COMMAND
        fi

        # 1) Unfold pseudodata with mW set to values of {0, 10, 20, ...}

        # for ((i=-100; i<=100; i+=10)); do
        for mass in "asimov" 100; do
            MABS=$mass
            if [ $mass -lt 0 ]; then
                UPDOWN="Down"
                MABS=$((-1*$MABS))
            elif [ $MABS -eq 0 ]; then
                UPDOWN=""
            else
                UPDOWN="Up"
            fi
            
            if [ "$mass" = "asimov" ]; then
                BIN=asimov
                echo "Perform unfolding for asimov data"
            else
                BIN=massShift${LABEL}${MABS}MeV$UPDOWN
                echo "Perform unfolding for $PSEUDO"
            fi

            FITRESULT=${COMBINE_ANALYSIS_OUTDIR}/fitresults_123456789_${BIN}.hdf5

            # 1) Unfold (pseudo)data
            if [ -e $FITRESULT ]; then
                echo "The file $FITRESULT exists, continue using it."
            elif [ "$BIN" = "asimov" ]; then
                COMMAND="cmssw-cc7 --command-to-run scripts/ci/setup_and_run_combine.sh $CMSSW_BASE $COMBINE_ANALYSIS_OUTDIR \
                    ${ANALYSIS}.hdf5 -t 1 --postfix $BIN --binByBinStat --correlateXsecStat"
            else
                PSEUDO="massWeight${LABEL}_massShift_${BIN}"
                COMMAND="cmssw-cc7 --command-to-run scripts/ci/setup_and_run_combine.sh $CMSSW_BASE $COMBINE_ANALYSIS_OUTDIR \
                    ${ANALYSIS}.hdf5 -p $PSEUDO --postfix $BIN --binByBinStat --correlateXsecStat"
                echo $COMMAND
                eval $COMMAND
            fi

            # 2)  Generate card with nominal theory model
            THEOMODEL_DIR=${COMBINE_OUTDIR}/${ANALYSIS}_${GEN_VARS_STR}${STATONLY}_${BIN}/
            THEOMODEL=${ANALYSIS}.hdf5
            THEOMODEL_PATH=${THEOMODEL_DIR}/${THEOMODEL}

            if [ -e $THEOMODEL_PATH ]; then
                echo "The file $THEOMODEL_PATH exists, continue using it."
            else
                COMMAND="./scripts/ci/run_with_singularity.sh scripts/ci/setup_and_run_python.sh scripts/combine/setupCombine.py \
                    -i $HISTMAKER_FILE -o $COMBINE_OUTDIR --hdf5 --fitvar $GEN_VARS $OPTTAU --postfix $BIN $OPTSTAT \
                    --fitresult $FITRESULT"
                echo $COMMAND
                eval $COMMAND        
            fi

            # 3) Fit theory model to unfolded pseudo data 
            FITRESULT_FINAL=${THEOMODEL_DIR}/fitresults_123456789.hdf5

            if [ -e $FITRESULT_FINAL ]; then
                echo "The file $FITRESULT_FINAL exists, continue using it."
            else
                COMMAND="cmssw-cc7 --command-to-run scripts/ci/setup_and_run_combine.sh $CMSSW_BASE $THEOMODEL_DIR $THEOMODEL --chisqFit --externalCovariance" 
                echo $COMMAND
                eval $COMMAND
            fi

            # 4) Compare observed mass pulls with pseudodata value

        done
    
    done

done
