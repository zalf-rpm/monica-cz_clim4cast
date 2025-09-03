#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

import asyncio
import capnp
import csv
from io import StringIO
import json
import numpy as np
import os
import sys
import zalfmas_fbp.run.components as c
import zalfmas_fbp.run.ports as p
import zalfmas_capnp_schemas

sys.path.append(os.path.dirname(zalfmas_capnp_schemas.__file__))
import fbp_capnp


async def run_component(port_infos_reader_sr: str, config: dict):
    ports = await p.PortConnector.create_from_port_infos_reader(
        port_infos_reader_sr, ins=["conf", "csv"], outs=["out"]
    )
    await p.update_config_from_port(config, ports["conf"])

    # get default country ids
    nuts3_ids = config["nuts3_ids"]

    while ports["csv"] and ports["out"]:
        try:
            msg = await ports["csv"].read()
            if msg.which() == "done":
                ports["csv"] = None
                continue
            csv_ip = msg.value.as_struct(fbp_capnp.IP)
            csv_str = csv_ip.content.as_text()

            observations = []
            with StringIO(csv_str) as file:
                dialect = csv.Sniffer().sniff(csv_str, delimiters=";,\t")
                file.seek(0)
                reader = csv.reader(file, dialect)
                next(reader, None)  # skip the header
                for row in reader:
                    id = int(row[24])
                    name = row[0]
                    for i in range(1, 24):
                        yield_t = float(row[i])
                        observations.append(
                            {
                                "id": id,
                                "nuts_id": name,
                                "year": 2000 + i - 1,
                                "value": np.nan
                                if yield_t < 0.0
                                else yield_t * 1000.0,  # t/ha -> kg/ha nan is -9999
                            }
                        )
            # order obs list by id to avoid mismatch between observation/evaluation lists
            observations.sort(key=lambda r: [r["id"], r["year"]])

            filtered_observations = list(
                filter(lambda d: d["id"] in nuts3_ids, observations)
            )
            if len(filtered_observations) == 0:
                continue

            param_set_id = "-".join([str(id) for id in nuts3_ids])
            out_ip = fbp_capnp.IP.new_message(
                attributes=[{"key": "param_set_id", "value": param_set_id}],
                content=json.dumps(filtered_observations)
            )
            await ports["out"].write(value=out_ip)

        except Exception as e:
            print(f"{os.path.basename(__file__)} Exception:", e)

    await ports.close_out_ports()
    print(f"{os.path.basename(__file__)}: process finished")


default_config = {
    "path_to_yield_data": "data/FAO_yield_data.csv",
    "nuts3_ids": [],
    "port:conf": "[TOML string] -> component configuration",
    "port:csv": "[string] -> csv string",
    "port:out": "",  # {country_id: {year: yield}} :string of json serialized mapping from country id to year to yield
}


def main():
    parser = c.create_default_fbp_component_args_parser(
        "Copy IP to all attached array out ports"
    )
    port_infos_reader_sr, config, args = c.handle_default_fpb_component_args(
        parser, default_config
    )
    asyncio.run(capnp.run(run_component(port_infos_reader_sr, config)))


if __name__ == "__main__":
    main()
