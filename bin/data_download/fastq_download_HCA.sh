#!/bin/bash

# download FASTQ files from human cell atlas data explorer from the paper: Cryopreservation and post-thaw characterization of dissociated human islet cells

curl --location --fail https://service.azul.data.humancellatlas.org/manifest/files/ksQwlKVkY3A0N6RjdXJsxBAN4eabSedfjrBmewOwVC2pxBBAvDEQ5bJcZJLSmdDwKyjtxCB6ZXBBBApHebC4--Kq8IkpgQ12QZ8UwlQUFmBdq-32Tg | curl --retry 15 --retry-delay 10 --config -