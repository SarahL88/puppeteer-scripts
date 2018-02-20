#!/usr/bin/env python
import uuid
import zipfile
import os
from Bio import SeqIO

import repository

OVERHANGSFILENAME = 'Overhangs.csv'
VECTORSFILENAME = 'Vectors.csv'
PLASMIDSFILENAME = 'Plasmids.csv'

PROTOCOL_FORMAT = 'edu-bu-synbiotools-format-moclo'
GENEFILE_FORMAT = 'gb'


def make_repo(archive, instanceid, authorid, date):
    repo = init_repository()
    set_default_format(repo);
    project = make_project(repo, instanceid, authorid, date);

    unzip_file(archive)
    archivefilepath = os.getcwd()

    subdirectoriesfolder, subdirectories = get_subdirectories(str(archive))
    overhangfiles, vectorfiles, plasmidfiles = get_filenames(subdirectories, subdirectoriesfolder)

    process_overhangs(repo, project, overhangfiles, instanceid, authorid, date);
    process_vectors(repo, project, vectorfiles, subdirectories, instanceid, authorid, date);
    process_plasmids(repo, project, plasmidfiles, subdirectories, instanceid, authorid, date);

    os.chdir(archivefilepath)

    return repo


def init_repository():
    repo = {}
    repo['collections'] = []
    repo['compositexrefs'] = []
    repo['cxref'] = []
    repo['families'] = []
    repo['features'] = []
    repo['ffxref'] = []
    repo['ffxrefpk'] = []
    repo['formats'] = []
    repo['gff3designs'] = []
    repo['nsa'] = []
    repo['nucseq'] = []
    repo['parts'] = []
    repo['plasmids'] = []
    repo['vectors'] = []
    return repo


def set_default_format(repo):
    format = {'idformat': PROTOCOL_FORMAT}
    repo['formats'].append(format)


def make_project(repo, instanceid, authorid, date):
    project = {}
    project['authorid'] = authorid
    project['datecreated'] = date
    project['description'] = 'Items described in project ' + instanceid
    project['name'] = 'project-' + instanceid
    project['idcollection'] = uuid.uuid4()
    repo['collections'].append(project)
    return project


def unzip_file(archive):
    with zipfile.ZipFile(archive, 'r') as archive:
        archive.extractall('HeadtoHead2.zip' + '-contents')


def get_subdirectories(archive):
    os.chdir(archive + '-contents/' + archive.strip('.zip'))
    subdirectoriesfolder = os.getcwd()
    subdirectories = os.listdir()
    if '.DS_Store' in subdirectories:
        subdirectories.remove('.DS_Store')
    return subdirectoriesfolder, subdirectories


def get_filenames(subdirectories, homefolder):
    overhangfiles = []
    vectorfiles = []
    plasmidfiles = []

    for subdir in subdirectories:
        os.chdir(subdir)
        files = os.listdir()
        overhangfiles.append([subdir + '/' + f for f in files if OVERHANGSFILENAME in f][0])
        vectorfiles.append([subdir + '/' + f for f in files if VECTORSFILENAME in f][0])
        plasmidfiles.append([subdir + '/' + f for f in files if PLASMIDSFILENAME in f][0])
        os.chdir(homefolder)

    return overhangfiles, vectorfiles, plasmidfiles


def process_overhangs(repo, project, overhangsfiles, instanceid, authorid, date):

    for overhangsfile in overhangsfiles:
        with open(overhangsfile) as file:
            inputlines = file.readlines()

        validate_input_file(inputlines[0], overhangsfile);

        if len(inputlines) <= 1:
            return

        collections = {}
        collections['authorid'] = authorid
        collections['datecreated'] = date
        collections['description'] ='Overhangs described in ' + instanceid
        collections['name'] = 'overhang-' + instanceid
        collectionid = uuid.uuid4()
        collections['idcollection'] = collectionid

        repo['collections'].append(collections)
        repository.add_object_to_collection(repo, collectionid, project['idcollection'], 'COLLECTION', authorid, date);

        for line in inputlines:
            if ',' not in line:
                continue;
            tokens = line.split(',')
            featurename = 'overhang-' + tokens[0]
            featuresequence = tokens[1].lower()
            featureid = repository.create_feature(repo, featurename, featuresequence, 'overhang', date);
            repository.add_object_to_collection(repo, collectionid, featureid, 'FEATURE', authorid, date);



def process_vectors(repo, project, vectorsfiles, directories, instanceid, authorid, date):

    for vectorsfile in vectorsfiles:
        with open(vectorsfile) as file:
            lines = file.readlines()
        validate_input_file(lines[0], vectorsfile);

        if len(lines) <= 1:
            continue

        vectors = {}
        vectors['authorid'] = authorid
        vectors['datecreated'] = date
        vectors['description'] = 'Vectors described in ' + instanceid
        collectionid = uuid.uuid4()
        vectors['idcollection'] = collectionid
        repo['collections'].append(vectors)
        repository.add_object_to_collection(repo, project['idcollection'],
                                         vectors['idcollection'], 'COLLECTION', authorid, date)

        lineno = 0
        for line in lines[1:]:
            tokens = line.split(',')
            tokens = [t.strip() for t in tokens if len(t.strip()) > 0]
            if len(tokens) < 5:
                raise ValueError('The Values.csv file does not have the required number of tokens on line ' + lineno)

            vectorfilename = tokens[0]
            vectorname = 'Vector-' + tokens[1]
            resistancename = 'Resistance-' + tokens[2]
            fiveprimeoverhangname = 'Overhang-' + tokens[3]
            threeprimeoverhangname = 'Overhang-' + tokens[4]

            if len(tokens) > 5:
                description = tokens[5]
            else:
                description = 'From ' + vectorfilename + ': ' + \
                              fiveprimeoverhangname + ', ' + \
                              threeprimeoverhangname + ', ' + resistancename;

            for directory in directories:
                files = os.listdir(directory)
                if vectorfilename in files:
                    vectorsequence = read_genbank_file(directory + '/' + vectorfilename)
                    break


            vector = {}
            vector['authorid'] = authorid
            vector['datecreated'] = date
            vector['description'] = description
            vector['name'] = vectorname
            vectorid = uuid.uuid4()
            vector['idvector'] = vectorid

            nucseq = {}
            nucseq['datecreated'] = date
            nucseq['idnucseq'] = vectorid
            nucseq['sequence'] = vectorsequence

            vector['nucseq'] = nucseq
            vector['iscircular'] = True

            repo['nucseq'].append(nucseq)
            repo['vectors'].append(vector)

            ft2 = {}
            foundfeature1 = False
            foundfeature2 = False

            overhangfeatures = repository.get_features_by_family_name(repo, 'overhang')

            for feature in overhangfeatures:

                if not foundfeature1 and feature['name'].upper() == fiveprimeoverhangname.upper():
                    position = repository.get_overhang_position_in_vector(vectorsequence, feature['nucseq']['sequence'])
                    repository.add_feature_to_nucseq(repo, vectorname, nucseq, feature, position, authorid, date)
                    foundfeature1 = True

                if not foundfeature2 and feature['name'].upper() == threeprimeoverhangname.upper():
                    position = repository.get_overhang_position_in_vector(vectorsequence, feature['nucseq']['sequence'])
                    repository.add_feature_to_nucseq(repo, vectorname, nucseq, feature, position, authorid, date)
                    ft2 = feature
                    foundfeature2 = True


            if not foundfeature1 or not foundfeature2:
                raise ValueError('The overhangs caused by vector ' + vectorname +
                                 ' were not defined in the overhangs manifest.')

            foundfeature = False

            for feature in repository.get_features_by_family_name(repo, 'resistance'):
                if feature['name'].upper() == resistancename.upper():
                    position = nucseq['sequence'].find(feature['nucseq']['sequence'])
                    repository.add_feature_to_nucseq(repo, vectorname, nucseq, feature, position, authorid, date)
                    foundfeature = True

            if not foundfeature:
                overhangpos = repository.get_overhang_position_in_vector(vectorsequence, ft2['nucseq']['sequence'])
                if overhangpos < 0:
                    raise ValueError('The overhang ' + ft2['name'] + ' could not be found in the vector ' + vectorname)
                startpos = overhangpos + len(ft2['nucseq']['sequence']) + 1
                resistancesequence = nucseq['sequence'][startpos: len(nucseq['sequence'])]
                f = repository.create_feature(repo, resistancename, resistancesequence, 'resistance', date)
                position = nucseq['sequence'].find(f['nucseq']['sequence'])
                repository.add_feature_to_nucseq(repo, vectorname, nucseq, f, position, authorid, date)
                repository.add_object_to_collection(repo, collectionid, f['idfeature'], 'FEATURE', authorid, date)

            repository.add_object_to_collection(repo, collectionid, vectorid, 'VECTOR', authorid, date)
            lineno += 1


def process_plasmids(repo, project, plasmidsfiles, directories, instanceid, authorid, date):

    for plasmidsfile in plasmidsfiles:
        with open(plasmidsfile) as file:
            lines = file.readlines()
        validate_input_file(lines[0], plasmidsfile);
        if len(lines) <= 1:
            continue

        plasmids = {}
        plasmids['authorid'] = authorid
        plasmids['datecreated'] = date
        plasmids['description'] = 'Parts described in ' + instanceid
        plasmids['name'] = 'part-' + instanceid
        collectionid = uuid.uuid4()
        plasmids['idcollection'] = collectionid
        repo['collections'].append(plasmids)

        # Create dict of 'family name : collection dict'
        allcollections = create_collections_by_family(repo, project, instanceid, authorid, date)

        lineno = 1
        for line in lines[1:]:
            tokens = line.split(',')
            tokens = [t.strip() for t in tokens if len(t.strip()) > 0]
            if len(tokens) < 4:
                raise ValueError('The plasmids.csv file does not comply with the format.'
                                 'Line ', lineno, ' has ', len(tokens),
                                 ' tokens, but should have at least 4.')

            plasmidfilename = tokens[0]
            partfamily = tokens[1].lower()
            partname = 'Part-' + tokens[2]
            vectorname = 'Vector-' + tokens[3]

            if len(tokens) > 4:
                description = tokens[4]
            else:
                description = 'From ' + plasmidfilename + ' part ' + partname + ' vector ' \
                              + vectorname + ' of type ' + partfamily + '.'

            for directory in directories:
                files = os.listdir(directory)
                if plasmidfilename in files:
                    plasmidsequence = read_genbank_file(directory + '/' + plasmidfilename)
                    break

            if not plasmidsequence:
                raise ValueError('Could not retrieve plasmid sequence.')

            partsequence = get_part_sequence(repo, plasmidsequence, vectorname)

            part = repository.persist_part(repo, partname, partsequence, description, True, authorid, date)
            persist_part_overhang_annotations(repo, vectorname, part, authorid, date)
            persist_part_feature(repo, vectorname, part, partfamily, authorid, date)

            repository.add_object_to_collection(repo, collectionid, part['idpart'], 'PART', authorid, date)

            # TODO - I add new families, instead of raising an error.  See Java line 702+
            if partfamily.lower() not in allcollections:
                add_new_family_to_all_collections(repo, partfamily, allcollections, project, instanceid, authorid, date)

            repository.add_object_to_collection(repo,
                                            allcollections[partfamily]['idcollection'], part['idpart'],
                                            'PART', authorid, date)

            vector = [v for v in repo['vectors'] if v['name'].lower() == vectorname.lower()][0]
            repository.persist_plasmid(repo, 'PLASMID-' + tokens[2], part, vector, authorid, date)

            lineno += 1



def create_collections_by_family(repo, project, instanceid, authorid, date):
    allcollections = {}
    families = repo['families']

    for fam in families:
        famname = fam['name']
        collection = {}
        collection['authorid'] = authorid
        collection['datecreated'] = date
        collection['description'] = 'Parts of ' + famname + ' family described in ' + instanceid
        collection['name'] = famname + '-part-' + instanceid
        collectionid = uuid.uuid4()
        collection['idcollection'] = collectionid
        repo['collections'].append(collection)
        allcollections[famname] = collection
        repository.add_object_to_collection(repo, project['idcollection'], collection['idcollection'], 'COLLECTION', authorid, date)

    return allcollections


def add_new_family_to_all_collections(repo, partfamily, allcollections, project, instanceid, authorid, date):
    family = {}
    family['name'] = partfamily
    family['idfamily'] = partfamily
    repo['families'].append(family)
    collection = {}
    collection['authorid'] = authorid
    collection['datecreated'] = date
    collection['description'] = 'Parts of ' + partfamily+ ' family described in ' + instanceid
    collection['name'] = partfamily + '-part-' + instanceid
    collectionid = uuid.uuid4()
    collection['idcollection'] = collectionid
    repo['collections'].append(collection)
    allcollections[partfamily] = collection
    repository.add_object_to_collection(repo, project['idcollection'], collection['idcollection'], 'COLLECTION', authorid, date)



def persist_part_overhang_annotations(repo, vectorname, part, authorid, date):
    nsa, fiveprimeoverhang, threeprimeoverhang = get_nucseq_annotations(repo, vectorname)

    repository.add_feature_to_nucseq(repo,
                                  fiveprimeoverhang['feature']['name'] + " in " + part['name'],
                                  part['nucseq'],
                                  fiveprimeoverhang['feature'],
                                  0,
                                  authorid,
                                  date)

    repository.add_feature_to_nucseq(repo,
                                  threeprimeoverhang['feature']['name'] + " in " + part['name'],
                                  part['nucseq'],
                                  threeprimeoverhang['feature'],
                                  len(part['nucseq']['sequence']) - len(threeprimeoverhang['feature']['nucseq']['sequence']),
                                  authorid,
                                  date)



def persist_part_feature(repo, vectorname, part, familyname, authorid, date):

    nsa, fiveprimeoverhang, threeprimeoverhang = get_nucseq_annotations(repo, vectorname)

    partseq = part['nucseq']['sequence'].strip().lower()
    start = len(fiveprimeoverhang['feature']['nucseq']['sequence'].strip())
    end = len(partseq) - len(threeprimeoverhang['feature']['nucseq']['sequence'].strip())
    partfeatureseq = partseq[start:end]

    partfeature = repository.create_feature(repo, 'Feature-' + part['name'], partfeatureseq, familyname, date)

    repository.add_feature_to_nucseq(repo,
                                  'Feature-' + part['name'],
                                  part['nucseq'],
                                  partfeature,
                                  len(fiveprimeoverhang['feature']['nucseq']['sequence']),
                                  authorid,
                                  date)


def get_part_sequence(repo, plasmidsequence, vectorname):
    nsa, fiveprimeoverhang, threeprimeoverhang = get_nucseq_annotations(repo, vectorname)

    fpopos = repository.get_overhang_position_in_plasmid(plasmidsequence,
                                                     fiveprimeoverhang['feature']['nucseq']['sequence'],
                                                     'FIVE_PRIME')

    tpopos = repository.get_overhang_position_in_plasmid(plasmidsequence,
                                                     threeprimeoverhang['feature']['nucseq']['sequence'],
                                                     'THREE_PRIME')

    partstart = fpopos
    partend = tpopos + len(threeprimeoverhang['feature']['nucseq']['sequence'].strip())-1
    partsequence = plasmidsequence[partstart:(partend+1)].strip()

    return partsequence.lower()



def get_nucseq_annotations(repo, vectorname):
    '''
    :return: nsa                nucseq annotation list   nsas associated with param vector
             fiveprimeoverhang  nucseq annotation
             threeprimeoverhang nucseq annotation
    '''
    nsa = get_vector_overhang_annotations(repo, vectorname)

    fiveprimeoverhang = nsa[0]
    threeprimeoverhang = nsa[1]

    if fiveprimeoverhang['startx'] > threeprimeoverhang['startx']:
        threeprimeoverhang = nsa[0]
        fiveprimeoverhang = nsa[1]

    return nsa, fiveprimeoverhang, threeprimeoverhang


def get_vector_overhang_annotations(repo, vectorname):

    vector = [v for v in repo['vectors'] if v['name'].lower() == vectorname.lower()]
    if len(vector) <= 0:
        raise ValueError('There is an invalid vector (' + vectorname + ') in the plasmids file.')
    vector = vector[0]

    nucseqannotations = repository.get_annotations_by_family(repo, vector['nucseq'], 'overhang')
    if len(nucseqannotations) != 2:
        raise ValueError('Unexpected number of overhangs found in ' + vectorname + ' num is ', len(nucseqannotations))

    return nucseqannotations


def read_genbank_file(filename):
   dnastring = ''
   for dna in SeqIO.parse(filename, GENEFILE_FORMAT):
       dnastring += str(dna.seq)

   return dnastring


def validate_input_file(firstline, filename):
    tokens = firstline.split(',')
    tokens = [t.lower() for t in tokens]
    valid = True

    if 'overhang' in filename:
        if (not tokens[0] == 'overhangname' or
                not tokens[1] == 'overhangsequence'):
            valid = False

    elif 'vector' in filename:
        if (not tokens[0] == 'filename' or
                not tokens[1] == 'vectorname' or
                not tokens[2] == 'resistance' or
                not tokens[3] == 'fiveprimeoverhang' or
                not tokens[4] == 'threeprimeoverhang' or
                not tokens[5] == 'description'):
            valid = False

    elif 'plasmid' in filename:
        if (not tokens[0] == 'filename' or
                not tokens[1] == 'partfamily' or
                not tokens[2] == 'partname' or
                not tokens[3] == 'vectorname' or
                not tokens[4] == 'description'):
            valid = False

    if not valid:
        raise ValueError('File ' + filename + ' is not in the correct format.')
