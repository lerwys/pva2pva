/*************************************************************************\
* Copyright (c) 2020 Michael Davidsaver
* EPICS BASE is distributed subject to a Software License Agreement found
* in file LICENSE that is included with this distribution.
\*************************************************************************/

#ifndef USE_TYPED_RSET
#  define USE_TYPED_RSET
#endif
#ifndef USE_TYPED_DSET
#  define USE_TYPED_DSET
#endif

#include <dbAccess.h>
#include <recSup.h>
#include <recGbl.h>
#include <dbEvent.h>
#include <errlog.h>
#include <alarm.h>

// include by ndaroiRecord.h
#include "epicsTypes.h"
#include "link.h"
#include "epicsMutex.h"
#include "ellLib.h"
#include "devSup.h"
#include "epicsTime.h"

#include "helper.h"

#include <pv/pvData.h>
#include <pv/standardField.h>
#include <pv/bitSet.h>
#include <pv/qsrv.h>

#define GEN_SIZE_OFFSET
#include <ndaroiRecord.h>
#undef  GEN_SIZE_OFFSET
#include <epicsExport.h>

#include <map>
#include <vector>
#include <utility>

namespace pvd = epics::pvData;

template <typename fldType>
struct fldCached {
    fldType *fld;   // record field address
    DBLINK *infld;  // record input link field address
    fldType val;    // cached record value for detecting changed

    fldCached(fldType *fld, DBLINK *infld, fldType val) : fld(fld), infld(infld), val(val) {};
    fldCached(fldType *fld, DBLINK *infld) : fld(fld), infld(infld), val(0) {};
};

struct ndaroiPvt {
    typedef std::map<std::string, std::tr1::shared_ptr<struct fldCached<epicsUInt32>> >
        mapType;
    typedef std::vector<mapType> vecType;
    // for each dimension:
    // * map of NTNTArray field name: record field offset
    vecType origDims;

    ndaroiPvt(vecType origDims) : origDims(origDims) {};
    ndaroiPvt() {};
};

namespace  {

ELLLIST vfPVStructureList = ELLLIST_INIT;
VFieldTypeNode vfPVStructureNode[3];

template <typename valueType>
static void convert_dim(const pvd::shared_vector<const valueType>& in,
        const pvd::PVStructureArrayPtr& in_dims,
        pvd::shared_vector<valueType>& out,
        const pvd::PVStructureArrayPtr& out_dims,
        int dim)
{
    epicsInt32 inStep, outStep, inOffset;
    epicsInt32 outSize;
    int incDir = 1;
    size_t inc;

    inOffset = (out_dims->view()[dim])->getSubFieldT<pvd::PVInt>("offset")->get();
    outSize = (out_dims->view()[dim])->getSubFieldT<pvd::PVInt>("size")->get();

    inStep = 1;
    outStep = 1;
    for (int i = 0; i < dim; ++i) {
        inStep  *= (in_dims->view()[i])->getSubFieldT<pvd::PVInt>("size")->get();
        outStep *= (out_dims->view()[i])->getSubFieldT<pvd::PVInt>("size")->get();
    }

    inc = incDir * inStep;

    pvd::shared_vector<const valueType> in_offs(in);
    in_offs.slice(inOffset*inStep);
    pvd::shared_vector<valueType> out_offs(out);

    for (epicsInt32 i = 0; i < outSize; ++i){
        if (dim > 0) {
            convert_dim<valueType>(in_offs, in_dims, out_offs, out_dims, dim-1);
        }
        else {
            out_offs[0] = in_offs[0];
        }

        in_offs.slice(inc);
        out_offs.slice(outStep);
    }
}

template <typename arrayType>
static void convert_dims(struct ndaroiRecord *prec, const pvd::PVStructureArrayPtr& new_dims)
{
    typedef typename arrayType::value_type valueType;
    typedef typename arrayType::const_svector vectorType;

    // array value
    std::tr1::shared_ptr<arrayType> pvarr = prec->val->getSubFieldT<pvd::PVUnion>("value")
                                                     ->get<pvd::PVValueArray<valueType> >();
    vectorType arr = pvarr->view();

    // dimensions
    pvd::PVStructureArrayPtr dims = prec->val->getSubFieldT<pvd::PVStructureArray>("dimension");

    // alloc buffer
    epicsInt32 length = 1;
    for (size_t i = 0; i < new_dims->getLength(); ++i) {
        length *= (new_dims->view()[i])->getSubFieldT<pvd::PVInt>("size")->get();
    }
    pvd::shared_vector<valueType> buf(length);

    convert_dim<valueType>(arr, dims, buf, new_dims, dims->getLength()-1);

    // update value/dimensions
    pvarr->replace(pvd::freeze(buf));
    prec->chg.set(pvarr->getFieldOffset());

    dims->copy(*new_dims);
    prec->chg.set(dims->getFieldOffset());
}

// convert dimensions
static void convert_dims(struct ndaroiRecord *prec, const pvd::PVStructureArrayPtr& new_dims)
{
    pvd::PVScalarArrayPtr pvarr = prec->val->getSubFieldT<pvd::PVUnion>("value")
                                     ->get<pvd::PVScalarArray>();

    // union can be null
    if (!pvarr) {
        throw std::runtime_error("empty union");
    }

    pvd::ScalarType stype = pvarr->getScalarArray()
                                 ->getElementType();

    switch(stype) {
        case pvd::pvBoolean: convert_dims<pvd::PVBooleanArray>(prec, new_dims); break;
        case pvd::pvByte:    convert_dims<pvd::PVByteArray>   (prec, new_dims); break;
        case pvd::pvShort:   convert_dims<pvd::PVShortArray>  (prec, new_dims); break;
        case pvd::pvInt:     convert_dims<pvd::PVIntArray>    (prec, new_dims); break;
        case pvd::pvLong:    convert_dims<pvd::PVLongArray>   (prec, new_dims); break;
        case pvd::pvUByte:   convert_dims<pvd::PVUByteArray>  (prec, new_dims); break;
        case pvd::pvUShort:  convert_dims<pvd::PVUShortArray> (prec, new_dims); break;
        case pvd::pvUInt:    convert_dims<pvd::PVUIntArray>   (prec, new_dims); break;
        case pvd::pvULong:   convert_dims<pvd::PVULongArray>  (prec, new_dims); break;
        case pvd::pvFloat:   convert_dims<pvd::PVFloatArray>  (prec, new_dims); break;
        case pvd::pvDouble:  convert_dims<pvd::PVDoubleArray> (prec, new_dims); break;
        default:
            throw std::runtime_error("unsupported array type");
            break;
    }
}

// convert PVScalarArray to the specified dimensions
static void convert(struct ndaroiRecord *prec, const pvd::PVStructureArrayPtr& new_dims)
{
    convert_dims(prec, new_dims);
}

long initialize()
{
    vfPVStructureNode[0].vtype = &vfStructure;
    ellAdd(&vfPVStructureList, &vfPVStructureNode[0].node);
    vfPVStructureNode[1].vtype = &vfPVStructure;
    ellAdd(&vfPVStructureList, &vfPVStructureNode[1].node);
    vfPVStructureNode[2].vtype = &vfSharedVector;
    ellAdd(&vfPVStructureList, &vfPVStructureNode[2].node);

    return 0;
}

long init_record(struct dbCommon *pcommon, int pass)
{
    ndaroiRecord *prec = (ndaroiRecord*)pcommon;
    ndaroidset *pdset = (ndaroidset *)(prec->dset);

    if(pass==0) {
        new (&prec->val) pvd::PVStructurePtr();
        prec->rpvt = new ndaroiPvt;
        new (&prec->chg) pvd::BitSet();
        new (&prec->vld) pvd::BitSet();
        new (&prec->ptyp) pvd::StructureConstPtr();

        if (!pdset) {
            recGblSetSevrMsg(prec, READ_ALARM, INVALID_ALARM, "no DSET");
            recGblRecordError(S_dev_noDSET, prec, "ndaroi: no DSET");
            return S_dev_noDSET;
        }

        // initialize dimension fields
        ndaroiPvt::mapType xdim, ydim, zdim;
        std::tr1::shared_ptr<struct fldCached<epicsUInt32>> xsize(new struct fldCached<epicsUInt32>(&prec->xsize, &prec->inxsize, 1));
        std::tr1::shared_ptr<struct fldCached<epicsUInt32>> xmin(new struct fldCached<epicsUInt32>(&prec->xmin, &prec->inxmin, 0));
        std::tr1::shared_ptr<struct fldCached<epicsUInt32>> ysize(new struct fldCached<epicsUInt32>(&prec->ysize, &prec->inysize, 1));
        std::tr1::shared_ptr<struct fldCached<epicsUInt32>> ymin(new struct fldCached<epicsUInt32>(&prec->ymin, &prec->inymin, 0));
        std::tr1::shared_ptr<struct fldCached<epicsUInt32>> zsize(new struct fldCached<epicsUInt32>(&prec->zsize, &prec->inzsize, 1));
        std::tr1::shared_ptr<struct fldCached<epicsUInt32>> zmin(new struct fldCached<epicsUInt32>(&prec->zmin, &prec->inzmin, 0));
        xdim["size"] = xsize;
        xdim["offset"] = xmin;
        ydim["size"] = ysize;
        ydim["offset"] = ymin;
        zdim["size"] = zsize;
        zdim["offset"] = zmin;
        prec->rpvt->origDims.push_back(xdim);
        prec->rpvt->origDims.push_back(ydim);
        prec->rpvt->origDims.push_back(zdim);

    } else { // pass==1
        if(pdset->common.init_record) {
            long ret = pdset->common.init_record(pcommon);
            if(ret)
                return ret;
        }

        if(!prec->val) {
            recGblRecordError(S_dev_noDSET, prec, "ndaroi: init_record must set VAL to a valid PVStructure");
            return S_db_badDbrtype;
        }

        prec->ptyp = prec->val->getStructure();
    }
    return 0;
}

long readValue(ndaroiRecord *prec, ndaroidset *pdset)
{
    epicsInt32 oldID = prec->id;
    long status = 0;

    /* fetch ROI dimensions */
    FOREACH(ndaroiPvt::vecType::iterator, it, end, prec->rpvt->origDims) {
        FOREACH(ndaroiPvt::mapType::iterator, it2, end2, *it) {

            // fetch link if not constant
            if(!status && !dbLinkIsConstant(it2->second->infld))
                status = dbGetLink(it2->second->infld, DBF_ULONG, it2->second->fld, 0, 0);

            // "size" cannot be 0, fix it
            if (!it2->first.compare("size") && *it2->second->fld == 0)
                *it2->second->fld = 1;

            // store cached value for detecting change
            it2->second->val = *it2->second->fld;
        }
    }

    if(status)
        return status;

    prec->chg.clear();

    // read pvstructure
    status = pdset->read_pvs(prec);

    // type change
    if(!status && prec->val && prec->val->getStructure() != prec->ptyp) {
        prec->ptyp = prec->val->getStructure();
    }

    try {
        // extract dims from new ROI dimensions
        pvd::PVStructureArrayPtr dims = prec->val->getSubFieldT<pvd::PVStructureArray>("dimension");
        pvd::StructureArrayConstPtr dims_type = dims->getStructureArray();
        pvd::PVStructureArrayPtr new_dims = dims_type->build();
        {
            pvd::StructureConstPtr elemType = dims_type->getStructure();
            pvd::PVStructureArray::svector arr;

            FOREACH(ndaroiPvt::vecType::iterator, it, end, prec->rpvt->origDims) {
                // build new dimension
                std::tr1::shared_ptr<pvd::PVStructure> dim = elemType->build();

                FOREACH(ndaroiPvt::mapType::iterator, it2, end2, *it)
                    dim->getSubFieldT<pvd::PVInt>(it2->first)->put(*it2->second->fld);

                arr.push_back(dim);
            }

            new_dims->replace(pvd::freeze(arr));
        }

        // convert to new dimensions
        convert(prec, new_dims);

            /* Post events if dimensions changed */
            FOREACH(ndaroiPvt::vecType::iterator, it, end, prec->rpvt->origDims) {
                FOREACH(ndaroiPvt::mapType::iterator, it2, end2, *it) {
                    if (it2->second->val != *it2->second->fld)
                        db_post_events(prec, it2->second->fld, DBE_VALUE|DBE_ARCHIVE);
                }
        }

        if(oldID != prec->id) {
            pvd::PVIntPtr id(prec->val->getSubFieldT<pvd::PVInt>("uniqueId"));
            id->put(prec->id);
            prec->chg.set(id->getFieldOffset());
        }

        prec->vld |= prec->chg;
    } catch(std::exception& e) {
        errlogPrintf("%s: readValue: %s\n", prec->name, e.what());
        recGblSetSevrMsg(prec, READ_ALARM, INVALID_ALARM, "readValue: %s", e.what());
        return -1;
    }

    return status;
}

void monitor(ndaroiRecord *prec)
{
    int monitor_mask = recGblResetAlarms(prec);

    if(monitor_mask & DBE_ALARM) {
        pvd::PVScalarPtr fld(prec->val->getSubFieldT<pvd::PVScalar>("alarm.severity"));
        fld->putFrom<pvd::uint16>(prec->sevr);
        prec->chg.set(fld->getFieldOffset());

        //TODO: map status properly
        fld = prec->val->getSubFieldT<pvd::PVScalar>("alarm.status");
        fld->putFrom<pvd::uint16>(prec->stat ? 1 : 0);
        prec->chg.set(fld->getFieldOffset());

        fld = prec->val->getSubFieldT<pvd::PVScalar>("alarm.message");
        fld->putFrom(std::string(prec->amsg));
        prec->chg.set(fld->getFieldOffset());
    }

    {
        pvd::PVScalarPtr fld(prec->val->getSubFieldT<pvd::PVScalar>("timeStamp.secondsPastEpoch"));
        fld->putFrom<pvd::uint32>(prec->time.secPastEpoch + POSIX_TIME_AT_EPICS_EPOCH);
        prec->chg.set(fld->getFieldOffset());

        fld = prec->val->getSubFieldT<pvd::PVScalar>("timeStamp.nanoseconds");
        fld->putFrom<pvd::uint32>(prec->time.nsec);
        prec->chg.set(fld->getFieldOffset());

        fld = prec->val->getSubFieldT<pvd::PVScalar>("timeStamp.userTag");
        fld->putFrom<pvd::uint64>(prec->utag);
        prec->chg.set(fld->getFieldOffset());

        monitor_mask |= DBE_VALUE|DBE_ARCHIVE;
    }

    db_post_events(prec, &prec->val, monitor_mask);
}

long process(struct dbCommon *pcommon)
{
    ndaroiRecord *prec = (ndaroiRecord*)pcommon;
    ndaroidset *pdset = (ndaroidset *)(prec->dset);
    unsigned char    pact=prec->pact;
    long status;

    if( (pdset == NULL) || (pdset->read_pvs == NULL) ) {
        prec->pact = TRUE;
        recGblSetSevrMsg(prec, READ_ALARM, INVALID_ALARM, "no read_pvs");
        errlogPrintf("%s: process: %s\n", prec->name, "read_pvs");
        return S_dev_missingSup;
    } else if(!prec->ptyp) {
        prec->pact = TRUE;
        recGblSetSevrMsg(prec, READ_ALARM, INVALID_ALARM, "no PTYP");
        return S_dev_NoInit;
    }

    status = readValue(prec, pdset); /* read the new value */
    /* check if device support set pact */
    if (!pact && prec->pact)
        return(0);

    prec->pact = TRUE;
    recGblGetTimeStamp(prec);

    try {
        /* check event list */
        monitor(prec);
    }
    catch (std::exception& e) {
        errlogPrintf("%s: process: %s\n", prec->name, e.what());
        recGblSetSevrMsg(prec, READ_ALARM, INVALID_ALARM, "process: %s", e.what());
        return S_db_badField;
    }

    /* process the forward scan link record */
    recGblFwdLink(prec);

    prec->pact=FALSE;
    return status;
}

long cvt_dbaddr(DBADDR *paddr)
{
    ndaroiRecord *prec = (ndaroiRecord*)paddr->precord;

    // we don't insist devsup allocate the required NELM
    paddr->ro = 1;
    // arbitrary limit
    paddr->no_elements = 1;
    // no access
    paddr->field_type = DBF_NOACCESS;

    // we provide vfield access for VAL and PTYP
    paddr->vfields = &vfPVStructureList;

    return 0;
}

long get_array_info(DBADDR *paddr, long *no_elements, long *offset)
{
    return S_db_badField;
}

long put_array_info(DBADDR *paddr, long nNew)
{
    return S_db_noMod;
}

long get_vfield(struct dbAddr *paddr, struct VField *p)
{
    ndaroiRecord *prec = (ndaroiRecord*)paddr->precord;

    try {
        if(p->vtype==&vfPVStructure) {
            VSharedPVStructure *pstr = (VSharedPVStructure*)p;
            if(dbGetFieldIndex(paddr)==ndaroiRecordVAL) {
                if(!*pstr->value)
                    return S_db_notInit;
                (*pstr->value)->copy(*prec->val);
                *pstr->changed = prec->vld;
                return 0;
            }

        } else if(p->vtype==&vfStructure) {
            VSharedStructure *pstr = (VSharedStructure*)p;
            if(dbGetFieldIndex(paddr)==ndaroiRecordVAL) {
                *pstr->value = prec->ptyp;
                return 0;
            }
        } else if(p->vtype==&vfSharedVector) {
            VSharedVector *pstr = (VSharedVector*)p;
            if(dbGetFieldIndex(paddr)==ndaroiRecordVAL) {
                pvd::PVScalarArrayPtr pvarr = prec->val->getSubFieldT<pvd::PVUnion>("value")
                                                       ->get<pvd::PVScalarArray>();
                if (!pvarr) {
                    throw std::runtime_error("empty \"value\" union");
                }
                pvarr->getAs<void>(*pstr->value);
                return 0;
            }
        }
    }
    catch (std::exception& e) {
        errlogPrintf("%s: get_vfield: %s\n", prec->name, e.what());
        return S_db_badChoice;
    }

    return S_db_badChoice;
}

long put_vfield(struct dbAddr *paddr, const struct VField *p)
{
    ndaroiRecord *prec = (ndaroiRecord*)paddr->precord;

    if(p->vtype==&vfPVStructure) {
        const VSharedPVStructure *pstr = (const VSharedPVStructure*)p;
        if(dbGetFieldIndex(paddr)==ndaroiRecordVAL) {
            prec->val->copy(**pstr->value);
            prec->vld |= *pstr->changed;
            return 0;
        }
    }
    return S_db_badChoice;
}

#define report NULL
#define special NULL
#define get_value NULL
#define get_units NULL
#define get_precision NULL
#define get_enum_str NULL
#define get_enum_strs NULL
#define put_enum_str NULL
#define get_graphic_double NULL
#define get_control_double NULL
#define get_alarm_double NULL

rset ndaroiRSET={
    RSETNUMBER,
    report,
    initialize,
    init_record,
    process,
    special,
    get_value,
    cvt_dbaddr,
    get_array_info,
    put_array_info,
    get_units,
    get_precision,
    get_enum_str,
    get_enum_strs,
    put_enum_str,
    get_graphic_double,
    get_control_double,
    get_alarm_double,
    get_vfield,
    put_vfield,
};

long readInitialPVStructure(struct link *pinp, void *raw)
{
    ndaroiRecord *prec = (ndaroiRecord *) pinp->precord;
    pvd::StructureConstPtr type;

    VSharedStructure ival;
    ival.vtype = &vfStructure;
    ival.value = &type;

    long status = dbGetLink(pinp, DBR_VFIELD, &ival, 0, 0);

    if (status)
        return status;

    prec->val = pvd::getPVDataCreate()->createPVStructure(type);

    return 0;
}

long readLocked(struct link *pinp, void *raw)
{
    const bool* doload = static_cast<const bool*>(raw);
    ndaroiRecord *prec = (ndaroiRecord *) pinp->precord;

    VSharedPVStructure ival;
    ival.vtype = &vfPVStructure;
    ival.value = &prec->val;
    ival.changed = &prec->chg;

    long status = *doload ? dbLoadLink(pinp, DBR_VFIELD, &ival) : dbGetLink(pinp, DBR_VFIELD, &ival, 0, 0);

    if (status)
        return status;

    if (dbLinkIsConstant(&prec->tsel) &&
        prec->tse == epicsTimeEventDeviceTime)
        dbGetTimeStamp(pinp, &prec->time);

    return status;
}

long init_record(struct dbCommon *pcommon)
{
    ndaroiRecord *prec = (ndaroiRecord *)pcommon;
    bool doload = true;

    if (readLocked(&prec->inp, &doload)) {
        // try to get at least Structure from link
        readInitialPVStructure(&prec->inp, 0);
        prec->udf = FALSE;
    }

    return 0;
}

long read_pvs(ndaroiRecord* prec)
{
    long status = 0;
    bool doload = false;

    if(!status)
        status = dbLinkDoLocked(&prec->inp, readLocked, &doload);

    if (status == S_db_noLSET)
        status = readLocked(&prec->inp, &doload);

    if (!status && !dbLinkIsConstant(&prec->inp))
        prec->udf = FALSE;

    return status;
}

ndaroidset devNDAROISoft = {
    {5, NULL, NULL, &init_record, NULL},
    &read_pvs
};

} // namespace

extern "C" {
epicsExportAddress(rset,ndaroiRSET);
epicsExportAddress(dset,devNDAROISoft);
}
