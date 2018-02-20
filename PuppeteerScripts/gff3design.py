import uuid

import repository
import fastapart


def make_GFF3_design(repo, design, countdesigns, partsdict):
    # TODO rename gff3designs and gff3parts - these are Constellation design results

    gff3design = {}
    gff3design['name'] = 'Design_' + str(countdesigns)
    gff3design['gff3parts'] = []

    parts = design.split(',')
    for part in parts:
        gff3part = {}
        gff3part['name'] = part
        gff3part['sequence'] = partsdict[part]
        gff3design['gff3parts'].append(gff3part)

    repo['gff3designs'].append(gff3design)
    return gff3design


def persist(repo, gff3parts, gff3design, authorid, datecreated):
    # Get parts in repo with the same name as constellation result parts ('gff3parts')
    constituentparts = get_constituent_parts(repo, gff3parts);

    #TODO - assess utility of making FASTAS.
    # We only use FASTAS to make compositexrefs objects,
    # which cross reference repo parts and constellation results.
    sequence = compute_part_sequence(repo, constituentparts)
    fastas = fastapart.make_fasta(repo, datecreated)
    part = fastapart.save_part(fastas, gff3design['name'], authorid, gff3design['name'], sequence, datecreated)

    save_constituent_part_references(repo, constituentparts, part)
    create_moclo_constituent_part_features(repo, part, constituentparts, authorid, datecreated)
    make_plasmid(repo, constituentparts, part, authorid, datecreated)



def get_constituent_parts(repo, gff3parts):
    constituentparts = []
    for gff3 in gff3parts:
        for part in repo['parts']:
            if gff3['name'].lower() in part['name'].lower():
                constituentparts.append(part)

    if not constituentparts:
        raise ValueError('Missing part returned in a combinatorial design.')

    return constituentparts


def save_constituent_part_references(repo, constituentparts, gff3design):
    '''
    :param repo:
    :param constituentparts: parts in repo that match constellation output
    :param part:             reference to FASTA part built based on constellation output
    :return:                 list of parts - for each constellation result part,
                                parent is constellation part
                                child is FASTA part
    '''
    count = 0
    for cp in constituentparts:
        compositexref = {}
        compositexref['childpart'] = gff3design    # Fasta part
        compositexref['parentpart'] = cp
        compositexref['idcomposite'] = uuid.uuid4()
        compositexref['position'] = count
        repo['compositexrefs'].append(compositexref)
        count +=1


def compute_part_sequence(repo, constituentparts):
    sequence = ''
    prev_endoverhang = ''
    count = 0
    for part in constituentparts:
        if 'moclo' in part['format']['idformat'].lower():

            # Get start and end overhangs for sequence
            overhangs = repository.get_annotations_by_family(repo, part['nucseq'], 'overhang')
            startoverhang = overhangs[0]
            endoverhang = overhangs[1]
            if overhangs[0]['startx'] > overhangs[1]['startx']:
                startoverhang = overhangs[1]
                endoverhang = overhangs[0]

            # Make sure the start overhang matches the previous sequence's end overhang
            if prev_endoverhang and startoverhang['feature']['nucseq']['sequence'] != \
                    prev_endoverhang['feature']['nucseq']['sequence']:
                raise ValueError("Consecutive composite parts don't have matching overhang values.")
            prev_endoverhang = endoverhang

            if startoverhang['startx'] > 0:
                raise ValueError('Part ' + part['name'] + ' does not begin with an overhang feature as expected '
                                                          'of a Modular Cloning part.')

            end_overhang_len = len(endoverhang['feature']['nucseq']['sequence'].strip())

            # Add start overhang to sequence
            sequence += startoverhang['feature']['nucseq']['sequence']\
                            .strip() + part['nucseq']['sequence'][:-end_overhang_len].strip()

            count += 1

    # Add end overhang to final sequence
    sequence += endoverhang['feature']['nucseq']['sequence'].strip()

    return sequence


def create_moclo_constituent_part_features(repo, pt, constituentparts, authorid, datecreated):
    '''
    Adds overhangs for constituent parts to repo.
    '''

    currfeaturestart = 0

    count = 0
    for cp in constituentparts:
        fiveprimeoverhang = repository.get_moclo_overhang_annotation(repo, cp['nucseq'], 'FIVE_PRIME')

        # For the first part, find the fiveprimeoverhang, add it to the repo
        if cp == constituentparts[0]:
            featurestart = pt['nucseq']['sequence'].strip()\
                .index(fiveprimeoverhang['feature']['nucseq']['sequence'].strip(), currfeaturestart)

            if featurestart < 0:
                raise ValueError ('Unable to find overhang annotation in the composite part.')
            if featurestart != currfeaturestart:
                raise ValueError("5' overhang found in composite part does not follow "
                                 "Moclo overhang rule: " + str(featurestart) + ' ' + str(currfeaturestart))
            repository.add_feature_to_nucseq(repo, pt['name'], pt['nucseq'],
                                          fiveprimeoverhang['feature'], featurestart, authorid, datecreated)
            currfeaturestart = featurestart + len(fiveprimeoverhang['feature']['nucseq']['sequence'].strip())


        # For each overhang annotation, find the start position in the sequence, and add it to the repo
        nsa = get_nonoverhang_annotations(repo, cp['nucseq'])
        for n in nsa:
            featurestart = pt['nucseq']['sequence'].strip().index(n['feature']['nucseq']['sequence']
                                                                  .strip(), currfeaturestart)

            if featurestart < 0:
                raise ValueError('Could not find ' + n['feature']['name'] + ' in ' + pt['name'] + '.')
            currfeaturestart = featurestart + len(n['feature']['nucseq']['sequence'].strip())
            repository.add_feature_to_nucseq(repo, pt['name'], pt['nucseq'],
                                          n['feature'], featurestart, authorid, datecreated)


        # Get the threeprimeoverhang anno, find where it starts, add it to the repo
        threeprimeoverhang = repository.get_moclo_overhang_annotation(repo, cp['nucseq'], 'THREE_PRIME')
        featurestart = pt['nucseq']['sequence'].strip()\
            .index(threeprimeoverhang['feature']['nucseq']['sequence'].strip(), currfeaturestart)

        if featurestart < 0:
            raise ValueError('Unable to find overhang annotation in the composite part.')
        if featurestart != currfeaturestart:
            raise ValueError("3' overhang found in composite part that does not follow Moclo rules: ", featurestart,
                             " ", currfeaturestart)
        repository.add_feature_to_nucseq(repo, pt['name'], pt['nucseq'], threeprimeoverhang['feature'],
                                      featurestart, authorid, datecreated)

        currfeaturestart = featurestart + len(fiveprimeoverhang['feature']['nucseq']['sequence'].strip())
        count += 1


def get_nonoverhang_annotations(repo, nucseq):
    allnucseqannos = [nsa for nsa in repo['nsa'] if nsa['nucseq']['idnucseq'] == nucseq['idnucseq']]

    nonoverhangannos = []
    for nsa in allnucseqannos:
        feature = nsa['feature']
        ffx = [ffx for ffx in repo['ffxref'] if ffx['feature']['idfeature'] == feature['idfeature']][0]
        fam = ffx['family']
        if 'overhang' not in fam['name'].lower():
            nonoverhangannos.append(nsa)
    return nonoverhangannos


def make_plasmid(repo, constituentparts, part, authorid, datecreated):
    first = constituentparts[0]
    last = constituentparts[len(constituentparts) - 1]

    nsafirst = repository.get_annotations_by_family(repo, first['nucseq'], 'overhang')
    nsalast = repository.get_annotations_by_family(repo, last['nucseq'], 'overhang')

    if nsafirst[0]['startx'] < nsafirst[1]['startx']:
        fpo = nsafirst[0]['feature']
    else:
        fpo = nsafirst[1]['feature']

    if nsalast[0]['startx'] < nsalast[1]['startx']:
        tpo = nsalast[0]['feature']
    else:
        tpo = nsalast[0]['feature']

    allresistances = repository.get_features_by_family_name(repo, 'resistance')

    resistancesforeachpart = get_resistance_by_composite_part(repo, constituentparts)
    availresistances = find_satisfiable_resistance(0, resistancesforeachpart, allresistances)

    if availresistances:
        eligiblevectors = repository.find_vectors_by_overhang_resistance(repo, fpo, tpo, availresistances)
        if eligiblevectors:
            return repository.persist_plasmid(repo, 'PLASMID-' +
                                             part['name'], part, eligiblevectors[0], authorid, datecreated)
    return {}


def get_resistance_by_composite_part(repo, constituentparts):
    resistancesforeachpart = []
    for pt in constituentparts:
        possiblevectors = repository.get_vectors_by_part(repo, pt)
        res = []
        for vt in possiblevectors:
            nsa = repository.get_annotations_by_family(repo, vt['nucseq'], 'resistance')
            if len(nsa) > 1:
                print('Warning: Multi resistance vector breaks assumption.')
            if nsa:
                res.append(nsa[0]['feature'])
        resistancesforeachpart.append(res)
    return resistancesforeachpart


def find_satisfiable_resistance(partnumber, currresforeachpart, availresistances):
    if partnumber < len(currresforeachpart):
        tochoosefrom = currresforeachpart[partnumber]
        for res in tochoosefrom:
            chosen = res
            if res in availresistances:
                remove = res
                availresistances.remove(res)
                result = find_satisfiable_resistance(partnumber + 1, currresforeachpart, availresistances)
                if not result:
                    availresistances.append(res)
                else:
                    return availresistances
    return availresistances
