import csv


def get_bunch(filein):
    fill_nums = {}
    with open(filein, "r") as lumicsv:
        reader = csv.reader(lumicsv)
        for row in reader:
            if row[0][0] == "#" or row[0][0] == "+":
                pass
            else:
                split_row = row[0].split("|")
                clean_row = [val for val in split_row if (val != "" and val != ",")]
                fill = clean_row[0].strip()
                if str(fill) not in fill_nums.keys():
                    nBunch = clean_row[-1].strip()
                    try:
                        fill_nums[str(fill)] = int(nBunch)
                    except:
                        pass
    return fill_nums


def add_bunch_col(filein, fileout, bunch_dict):
    with open(filein, "r") as lumicsv:
        with open(fileout, "w") as outfile:
            writer = csv.writer(outfile)
            reader = csv.reader(lumicsv)
            for row in reader:
                # print(row)
                if row[0][0:2] == "#D":
                    writer.writerow(row)
                    # pass
                elif row[0][0:2] == "#r":
                    # pass
                    writer.writerow(row + ["nBunches"])
                else:
                    try:
                        _, fill = row[0].split(":")
                        n_bunches = bunch_dict[fill]
                        writer.writerow(row + [str(n_bunches)])
                    except:
                        pass


nBunch_dict = get_bunch("/work/submit/jbenke/WRemnants/wremnants-data/data/byfill.csv")
print(nBunch_dict)
add_bunch_col(
    "/work/submit/jbenke/WRemnants/wremnants-data/data/bylsoutput.csv",
    "/work/submit/jbenke/WRemnants/wremnants-data/data/bylsoutput_nBunches.csv",
    nBunch_dict,
)
add_bunch_col(
    "/work/submit/jbenke/WRemnants/wremnants-data/data/bylsoutput_PCC.csv",
    "/work/submit/jbenke/WRemnants/wremnants-data/data/bylsoutput_nBunches_PCC.csv",
    nBunch_dict,
)
add_bunch_col(
    "/work/submit/jbenke/WRemnants/wremnants-data/data/bylsoutput_HFOC.csv",
    "/work/submit/jbenke/WRemnants/wremnants-data/data/bylsoutput_nBunches_HFOC.csv",
    nBunch_dict,
)
add_bunch_col(
    "/work/submit/jbenke/WRemnants/wremnants-data/data/bylsoutput_RAMSES.csv",
    "/work/submit/jbenke/WRemnants/wremnants-data/data/bylsoutput_nBunches_RAMSES.csv",
    nBunch_dict,
)


# def add_bunch_col(filein):
#     fill_nums = []
#     with open(filein, 'r') as lumicsv:
#         reader = csv.reader(lumicsv)
#         for row in reader:
#             if row[0][0]=="#":
#                 continue
#             else:
#                 _, fill = row[0].split(":")
#                 if int(fill) not in fill_nums:
#                     fill_nums.append(int(fill))
#     return fill_nums


# print(add_bunch_col('/work/submit/jbenke/WRemnants/wremnants-data/data/bylsoutput_nBunches.csv'))


# import subprocess
# import csv


# fill_cache = {}

# def get_number_of_bunches(fill_number):
#     if fill_number in fill_cache:
#         return fill_cache[fill_number]
#     cmd = [
#         "singularity", "exec",
#         "--bind", "$HOME,$HOME/.brilconda,/cvmfs",
#         "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cloud/brilws-docker:latest",
#         "brilcalc", "fill", "-f", str(fill_number), "--output-style", "csv"
#     ]

#     proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

#     stdout, stderr = proc.communicate()

#     print("brilcalc output:")
#     print(stdout)
#     if not stdout.strip():
#         print("Empty stdout from brilcalc!")
#         print("stderr:", stderr.strip())
#     lines = stdout.strip().split("\n")

#     # print("header_line:", header_line)
#     # print("type(header_line):", type(header_line))
#     header_line = None
#     data_line = None
#     for line in lines:
#         print(line)
#         if line.startswith("#"):
#             header_line = line.lstrip("#").strip()
#         elif header_line and not data_line:
#             data_line = line.strip()
#             break
#     print("header_line:", header_line)
#     print("type(header_line):", type(header_line))
#     headers = [h.strip() for h in header_line.split(",")]
#     values = [v.strip() for v in data_line.split(",")]
#     n_bunches = int(values[headers.index("nCollidingBunches")])

#     return n_bunches

# def add_bunch_col(filein, fileout):
#     with open(filein, 'r') as lumicsv:
#         with open(fileout, 'w') as outfile:
#             writer = csv.writer(outfile)
#             reader = csv.reader(lumicsv)
#             for row in reader:
#                 # print(row)
#                 if row[0][0:2]=="#D":
#                     writer.writerow(row)
#                     # continue
#                 elif row[0][0:2] =='#r':
#                     # row.append(",nBunches")
#                     writer.writerow(row + ['nBunches'])
#                 else:
#                     _, fill = row[0].split(":")
#                     n_bunches = get_number_of_bunches(fill)
#                     writer.writerow(row + [str(n_bunches)])


# add_bunch_col('bylsoutput_nBunches.csv', 'bylsoutput_nBunches_out.csv')
