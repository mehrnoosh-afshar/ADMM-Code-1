/* Copyright 2015 The MathWorks, Inc. */

#ifndef PTP_TYPES_H
#define PTP_TYPES_H

#include <tmwtypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <windows.h>
#include <stdint.h>
#include <math.h>

#pragma pack(push,1)

#define MAX_MESSAGE_LENTGH     256

typedef enum {
    PTP_NO_ERROR,
    PTP_DATATYPE_OVERFLOW,
    PTP_DATA_TYPE_UNDERFLOW,
    PTP_NULL_POINTER,
    PTP_UNKNOWN_PARAM,
    PTP_NO_TIMESTAMP,
    PTP_UNMATCHED_TIMESTAMP,
    PTP_BOARD_ERROR,
    PTP_READ_HW_TIME_ERR,
    PTP_UNSUPPORTED_OPTION,
    PTP_IRELEVANT_MSG_RX,
    PTP_INCOMPATIBLE_OPTIONS,
    PTP_WIN_FUNC_ERR,
    PTP_INVALID_ARG,
    PTP_INVALID_DATA,
    PTP_FAULTY_STATE,
    PTP_PROTOCOL_REQ_ERROR,
    PTP_DATASET_COMP_ERROR,
    PTP_BMC_ERROR,
    PTP_BOARD_NOT_SUP,
    PTP_UNKNOWN_ERROR
} PTP_ReturnType;

typedef enum {
    INITIALIZING = 0x01,
    FAULTY,
    DISABLED,
    LISTENING,
    PRE_MASTER,
    MASTER,
    PASSIVE,
    UNCALIBRATED,
    SLAVE
} PTP_PortState;

typedef enum {
    GET_MNGMT,
    SET_MNGMT,
    RESPONSE_MNGMT,
    COMMAND_MNGMT,
    ACKNOWLEDGE_MNGMT
} PTP_ManagementActionField;

typedef enum {
   
    NULL_MANAGEMENT,
    CLOCK_DESCRIPTION_MNGMT,
    USER_DESCRIPTION_MNGMT,
    SAVE_IN_NON_VOLATILE_STORAGE_MNGMT,
    RESET_NON_VOLATILE_STORAGE_MNGMT,
    INITIALIZE_MNGMT,
    FAULT_LOG_MNGMT,
    FAULT_LOG_RESET_MNGMT,
   
    DEFAULT_DATA_SET_MNGMT = 0x2000,
    CURRENT_DATA_SET_MNGMT,
    PARENT_DATA_SET_MNGMT,
    TIME_PROPERTIES_DATA_SET_MNGMT,
    PORT_DATA_SET_MNGMT,
    PRIORITY1_MNGMT,
    PRIORITY2_MNGMT,
    DOMAIN_MNGMT,
    SLAVE_ONLY_MNGMT,
    LOG_ANNOUNCE_INTERVAL_MNGMT,
    ANNOUCE_RECEIPT_TIMEOUT_MNGMT,
    LOG_SYNC_INTERVAL_MNGMT,
    VERSION_NUMBER_MNGMT,
    ENABLE_PORT_MNGMT,
    DISABLE_PORT_MNGMT,
    TIME_MNGMT,
    CLOCK_ACCURACY_MNGMT,
    UTC_PROPERTIES_MNGMT,
    TRACEABILITY_PROPERTIES_MNGMT,
    TIMESCALE_PROPERTIES_MNGMT,
    UNICAST_NEGOTIATION_ENABLE_MNGMT,
    PATH_TRACE_LIST_MNGMT,
    PATH_TRACE_ENABLE_MNGMT,
    GRANDMASTER_CLUSTER_TABLE_MNGMT,
    UNICAST_MASTER_TABLE_MNGMT,
    UNICAST_MASTER_MAX_TABLE_SIZE_MNGMT,
    ACCEPTABLE_MASTER_TABLE_MNGMT,
    ACCEPTABLE_MASTER_TABLE_ENABLED_MNGMT,
    ACCEPTABLE_MASTER_MAX_TABLE_SIZE_MNGMT,
    ALTERNATE_MASTER_MNGMT,
    ALTERNATE_TIME_OFFSET_ENABLE_MNGMT,
    ALTERNATE_TIME_OFFSET_NAME_MNGMT,
    ALTERNATE_TIME_OFFSET_MAX_KEY_MNGMT,
    ALTERNATE_TIME_OFFSET_PROPERTIES_MNGMT,
   
    TRANSPARENT_CLOCK_DEFAULT_DATA_SET_MNGMT = 0x4000,
    TRANSPARENT_CLOCK_PORT_DATA_SET_MNGMT,
    PRIMARY_DOMAIN_MNGMT,
   
    DELAY_MECHANISM_MNGMT = 0x6000,
    LOG_MEAN_PDELAY_REQ_INTERVAL_MNGMT
} PTP_ManagementId;

typedef enum {
    MANAGEMENT_TLV = 1,
    MANAGEMENT_ERROR_STATUS,
    ORGANIZATION_EXTENSION,
    REQUEST_UNICAST_TRANSMISSION,
    GRANT_UNICAST_TRANSMISSION,
    CANCEL_UNICAST_TRANSMISSION,
    PATH_TRACE,
    ALTERNATE_TIME_OFFSET_INDICATOR,
    AUTHENTICATION = 0x2000,
    AUTHENTICATION_CHALLENGE,
    SECURITY_ASSOCIATION_UPDATE,
    CUM_FREQ_SCALE_FACTOR_OFFSET
} PTP_TlvTypes;

typedef struct {
    uint8_T       deviceIndex;
    uint32_T      board_time_increment;
    uint8_T       time_source;
    uint8_T       time_scale;
    SYSTEMTIME    arbitrary_epoch;
    uint8_T       delay_mech;
    uint8_T       domain_number;
    uint8_T       priority1;
    uint8_T       priority2;
    uint8_T       clockClass;
    uint8_T       clockAccuracy;
    boolean_T     slave_only;
    real_T        sync_interval;
    real_T        announce_interval;
    real_T        min_pdelay_req_interval;
    uint32_T      announce_receipt_timeout;
    real_T        sample_time;
    uint32_T      ip_address;
    uint8_T       transport_protocol;
    uint16_T      utc_offset;
    int32_T       bus;
    int32_T       slot;
    int32_T       function;
} PTP_ConfigParams;

typedef struct {
    boolean_T                       isNewMessage;
    uint16_T                        sequenceNumber;
    PTP_ManagementId                managementId;
    PTP_ManagementActionField       action;
    PTP_TlvTypes                    tlvType;
    uint16_T                        tlvDataLength;
} PTP_ManagementMessageInfo;

#pragma pack(pop)

#endif