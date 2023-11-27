#!/usr/bin/env python

"""
String formatting to enter shadie.Model.map dictionary events
information into ordered SLiM3 code blocks.
"""

from typing import Dict, List, Union

INITIALIZE = """
initialize() {{

  // model type
  initializeSLiMModelType("nonWF");

  // config
  initializeRecombinationRate({recombination_rate});
  initializeMutationRate({mutation_rate});
  initializeTreeSeq(simplificationInterval={simplification_interval});

  // MutationType init
  {mutations}
  c({mutation_names}).haploidDominanceCoeff = 1.0;

  // ElementType init
  {elements}

  // Chromosome (GenomicElement init)
  {chromosome}

  // constants (Population sizes and others)
  {constants}

  // globals (metadata dictionary)
  {simglobals}

  // extra scripts (Optional)
  {scripts}
}}
"""
# --------------------------------------------

REPRODUCTION = """
// generates offspring
{comment}{idx}reproduction({population}) {{
    {scripts}
}}
"""

MUT_EFFECT = """
// adjusts mutationEffect calculation 
//(specific mutation in specific ind)
{comment}{idx}mutationEffect({mutation}) {{
    {scripts}
}}
"""

FITNESS = """
// adjusts fitnessEffect calculation 
//(fitness recalc for specific ind)
{comment}{idx}fitnessEffect({population}) {{
    {scripts}
}}
"""

SURVIVAL = """
// implements survival adjustments
{comment}{idx}survival({population}) {{
    {scripts}
}}
"""

FIRST = """
// executes before offspring are generated
{comment}
{idx}{time}first() {{
    {scripts}
}}
"""


EARLY = """
// executes after offspring are generated
{comment}
{idx}{time}early() {{
    {scripts}
}}
"""


LATE = """
// executes after selection occurs
{comment}
{idx}{time}late() {{
    {scripts}
}}
"""

CUSTOM = """{comment}{scripts}"""

EVENT_TO_FORMATTER = {
    "initialize": INITIALIZE,
    "first": FIRST,
    "early": EARLY,
    "late": LATE,
    'muteffect': MUT_EFFECT,
    'fitness': FITNESS,
    'survival': SURVIVAL,
    'custom': CUSTOM,
    'reproduction': REPRODUCTION,
    'shadie': CUSTOM,
}


def clean_scripts(scripts: Union[str, List[str]]):
    """
    Ensures scripts end with a semi-colon
    """
    if isinstance(scripts, list):
        scripts = "\n    ".join([i.strip(';') + ';' for i in scripts])
    else:
        scripts = scripts.strip()
        scripts = (
            scripts if scripts.endswith("}") else scripts.strip(";") + ";"
        )
    return scripts


def format_event_dicts_to_strings(event: Dict):
    """
    Performs string formatting on the .map dictionary of the
    Model object to write the arguments to SLiM string format
    """
    # cleanup formatting of some arguments
    if 'constants' in event:
        event['constants'] = "\n  ".join([
            f"defineConstant('{str(key).upper()}', {val});" for key, val
            in event['constants'].items()
        ])

    if 'simglobals' in event:
        for (key, val) in event['simglobals'].items():
            if isinstance(val, dict):
                dictlist = []
                for (i,j) in val.items():
                    dictlist.append(str(i))
                    dictlist.append(str(j))
                string = "Dictionary"+str(tuple(dictlist))
                event['simglobals'].update({key:str(string)})

        event['simglobals'] = "\n  ".join([
            f"defineGlobal('{str(key).upper()}', {val});" for key, val
            in event['simglobals'].items()
        ])

    if 'scripts' in event:
        event['scripts'] = clean_scripts(event['scripts'])

    if 'time' in event:
        event['time'] = f"{event['time']} " if event['time'] else ""

    if 'comment' in event:
        event['comment'] = (
            "// " + event['comment'].lstrip("//").strip() + "\n"
            if event['comment'] else "")

    if 'idx' in event:
        event['idx'] = f"{event['idx']} " if event['idx'] else ""

    if 'mutation' in event:
        event['mutation'] = event.get("mutation", "")

    if 'population' in event:
        event['population'] = (
            f"{event['population']} " if event['population'] else "")
    return event
