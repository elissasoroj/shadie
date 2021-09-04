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
  initializeTreeSeq(checkCoalescence=T);

  // MutationType init
  {mutations}

  // ElementType init
  {elements}

  // Chromosome (GenomicElement init)
  {chromosome}

  // constants (Population sizes and others)
  {constants}

  // extra scripts (Optional)
  {scripts}
}}
"""
# --------------------------------------------

REPRODUCTION = """
// generates offspring
{comment}reproduction({population}) {{
    {scripts}
}}
"""

FITNESS = """
// adjusts fitness calculation
{comment}{idx}fitness({mutation}) {{
    {scripts}
}}
"""

SURVIVAL = """
// implements survival adjustments
{idx}survival({population}) {{
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
    "early": EARLY,
    "late": LATE,
    'fitness': FITNESS,
    'survival': SURVIVAL,
    'custom': CUSTOM,
    'reproduction': REPRODUCTION,
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
            f"defineConstant('{key}', {val});" for key, val
            in event['constants'].items()
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
