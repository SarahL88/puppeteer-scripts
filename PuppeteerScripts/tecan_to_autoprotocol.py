#!/usr/bin/env python
import sys
import json
import collections
from autoprotocol.protocol import Protocol # pip install autoprotocol

# Output
OUTPUT_AUTOPROTOCOL_FILE = 'autoprotocol.json'

SOURCEPLATE = 'DNASourcePlate'
DESTPLATE = 'MoCloDestinationPlate'
WELL = '96 Well Microplate'

MAX_SOURCEPLATE_ROW = 6
MAX_DESTPLATE_ROW   = 7

# row/col -> well number
rc_to_wn = {}

# name -> well number
names_to_wells = {}

def convert_tecan(jsonRequest):

    tecan_directions = jsonRequest['tecanProgram']

    # Get text block (from 'aspirate' line to 'droptips' line)
    txtblock, rest_of_file = get_aspirate_to_droptips(tecan_directions)
    #print(txtblock)

    # Dict of well number : part
    wells_parts_dict = get_parts_dict(tecan_directions)

    # Dict of well number : line that identifies this well number
    aspirate_dict = {}
    aspirate_wells = []
    dispense_wells = []
    volume_list = []

    # Source/dest commands
    aspirate = []
    dispense = []

    # Source/dest well numbers
    aspirate_nums = []
    dispense_nums = []
    volumes = []

    # For each chunk of text, get aspirate/dispense lines
    while txtblock:

        for line in txtblock.splitlines():
            tokens = line.split(',')
            row = int(tokens[4].split('=')[-1])
            col = int(tokens[5].split('=')[-1])
            volume = tokens[7].split('=')[1].strip()

            volm = int(volume.split('.')[0])
            volume_list.append(volm)
            if volm > 2:
                if volm not in volumes:
                    volumes.append(volm)

            if 'aspirate' in tokens[0].lower():
                well_number = get_well_number(row, col, 'source')
                aspirate_command = 'A;' + SOURCEPLATE + ';;' + WELL + ';' + str(well_number) + ';;' + volume
                aspirate.append(aspirate_command)
                aspirate_dict[well_number] = aspirate_command
                aspirate_wells.append(well_number)
                if well_number not in aspirate_nums:
                    aspirate_nums.append(well_number)

            elif 'dispense' in tokens[0].lower():
                well_number = get_well_number(row, col, 'dest')
                dispense.append('D;' + DESTPLATE + ';;' + WELL + ';' + str(well_number) + ';;' + volume)
                dispense_wells.append(well_number)
                if well_number not in dispense_nums:
                    dispense_nums.append(well_number)

        txtblock, rest_of_file = get_aspirate_to_droptips(rest_of_file)

    instructions_dict = combine_tecan_instructions(wells_parts_dict, aspirate_wells, dispense_wells, volume_list)
    return generate_autoprotocol(instructions_dict)

def get_well_number(row, col, plate_type):
    if 'source' in plate_type:
        wn = row + (( col - 1) * MAX_SOURCEPLATE_ROW)

    elif 'dest' in plate_type:
        wn = row + ((col-1) * MAX_DESTPLATE_ROW)
    else:
        raise ValueError("Wrong plate type provided.  Plate type must be 'source' or 'dest'.")
    return wn


def get_aspirate_to_droptips(txt):
    # Get chunk of text from 'aspirate' line to 'droptips' line
    aspirateindex = txt.find('aspirate')
    if aspirateindex == -1:
        return '', ''
    txt = txt[aspirateindex:]
    droptips_index= txt.find('dropTips')
    return txt[:droptips_index], txt[(droptips_index-1):]


def get_parts_dict(tecan_directions):
    wells_parts_dict = {}  # well_number : part name
    tecan_directions = tecan_directions[:tecan_directions.find('aspirate')]

    well_numbers = []
    part_name = ''
    count_master_mix = 0
    for line in tecan_directions.splitlines():
        if 'Part-' in line:
            part = line.split('Part-')[1].split(' in well')[0]
            part_name = part.split('-')[0]
        elif 'Vector-' in line:
            part = line.split('Vector-')[1].split(' in well')[0]
            part_name = part.split('-')[0]
        elif 'Master Mix' in line:
            part_name = 'Master Mix'
        elif 'backbone' in line:
            part_name = 'backbone'
        if len(part_name) > 0 and 'in well' in line:
            well_num = line.split('in well ')[1].split(' of plate')[0]
            row = int(ord(well_num[0].lower()) - 96)
            col = int(well_num[1:])
            rowcol = str(row)+well_num[1:]
            if 'Master' in part_name:
                if count_master_mix > 0:
                    part_name = part_name + str(count_master_mix)
                rc_to_wn[rowcol] = get_well_number(row, col, 'source')
                count_master_mix += 1
            else:
                rc_to_wn[rowcol] = get_well_number(row, col, 'source')

            well_number = get_well_number(row, col, 'source')
            well_numbers.append(well_number)
            wells_parts_dict[well_number] = part_name
            names_to_wells[part_name] = well_number
            part_name = ''

    well_numbers.sort()
    return collections.OrderedDict(sorted(wells_parts_dict.items()))

def combine_tecan_instructions(parts_dict, aspirate_wells, dispense_wells, volumes):
    instructions_dict = {}
    instructions_dict['parts'] = parts_dict
    instructions_dict['source_wells'] = aspirate_wells
    instructions_dict['dest_wells'] = dispense_wells
    instructions_dict['volumes'] = volumes
    return instructions_dict

def generate_autoprotocol(instructions_dict):
    p = Protocol()
    part_ref_dict = {} # ref name: ref_object
    # add destination plate to ref dict
    moclo_dest_plate_name = 'MoCloDestinationPlate'
    moclo_dest_plate_ref = p.ref(moclo_dest_plate_name, cont_type="96-pcr", discard=True)
    part_ref_dict[moclo_dest_plate_name] = moclo_dest_plate_ref

    # create refs
    for num,part_name in instructions_dict['parts'].items():
        part_ref = p.ref(part_name, cont_type="96-pcr", discard=True)
        part_ref_dict[part_name] = part_ref

    # create instructions
    for x in range(0, len(instructions_dict['source_wells'])):
        from_well = instructions_dict['source_wells'][x]
        part_name = instructions_dict['parts'][from_well] # get part name in that well
        source_plate = part_ref_dict[part_name] # get ref object
        volume = instructions_dict['volumes'][x] # get volume
        volume_str = str(volume) + ':microliter'
        to_well = instructions_dict['dest_wells'][x] # get destination well num

        p.transfer(source_plate.wells_from(from_well,1), moclo_dest_plate_ref.wells_from(to_well,1),volume_str)

    return p.as_dict()


# LOCAL TESTING
# if __name__ == "__main__":
#     with open('pup_response.json') as json_data:
#         d = json.load(json_data)
#         protocol_dict = convert_tecan(d)
#
#         print(json.dumps(protocol_dict))
