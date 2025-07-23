"""
Define classes to represent entities in metadata
"""
import copy


class MappingMetadata(object):

    def __init__(
            self,
            src,
            dst,
            citation):

        if not isinstance(src, MetadataEntity):
            raise ValueError(
                'MappingMetadata.src must be a MetadataEntity; '
                f'you gave {src} which is a {type(src)}'
            )

        if not isinstance(dst, MetadataEntity):
            raise ValueError(
                'MappingMetadata.dst must be a MetadataEntity; '
                f'you gave {dst} which is a {type(dst)}'
            )

        if not isinstance(dst, type(src)):
            raise ValueError(
                'MappingMetadata.src and dst must be of '
                f'same types. You gave {type(src)} and {type(dst)}'
            )

        if not isinstance(citation, dict):
            raise ValueError(
                'MappingMetadata.citation must be of type dict; '
                f'you gave {type(citation)}'
            )

        if isinstance(src, Authority):
            axis = 'authority'
        elif isinstance(src, Species):
            axis = 'species'
        else:
            raise ValueError(
                f'Cannot infer axis from src of type {type(src)}'
            )

        self._axis = axis
        self._citation = copy.deepcopy(citation)

        if 'idx' in self._citation:
            self._citation.pop('idx')

        self._src = src
        self._dst = dst

    def serialize(self):
        return {
            'mapping': {
                'axis': self._axis,
                'from': self._src.serialize(),
                'to': self._dst.serialize()
            },
            'citation': self._citation
        }


class MetadataEntity(object):
    pass


class Authority(MetadataEntity):

    def __init__(
            self,
            name):

        valid = ('ENSEMBL', 'NCBI', 'symbol')
        if name not in valid:
            raise ValueError(
                f"authority = '{name}' is not a valid authority; "
                f"must be one of {valid}"
            )

        if not isinstance(name, str):
            raise ValueError(
                'Authority.name must be a str; you gave '
                f'{name} which is of type {type(name)}'
            )
        self._name = name

    def __repr__(self):
        return f"{self._name}"

    def serialize(self):
        return f'{self._name}'

    @property
    def name(self):
        return self._name


class Species(MetadataEntity):

    def __init__(
            self,
            name,
            taxon):

        if not isinstance(name, str):
            raise ValueError(
                'Species.name must be a str; you gave '
                f'{name} which is of type {type(name)}'
            )

        if not isinstance(taxon, int):
            raise ValueError(
                'Species.taxon must be a int; you gave '
                f'{taxon} which is of type {type(taxon)}'
            )

        self._name = name
        self._taxon = taxon

    def __repr__(self):
        return f"{self._name}:{self._taxon}"

    def serialize(self):
        return {
            'name': self._name,
            'taxon': self._taxon
        }

    @property
    def name(self):
        return self._name

    @property
    def taxon(self):
        return self._taxon
