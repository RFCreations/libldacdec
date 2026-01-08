#include <assert.h>
#include <string.h>

#include "ldacdec.h"
#include "log.h"
#include "utility.h"
#include "bit_reader.h"
#include "imdct.h"
#include "huffCodes.h"
#include "spectrum.h"
#include "bit_allocation.h"

#define LDAC_SYNCWORDBITS   (8)
#define LDAC_SYNCWORD       (0xAA)
/** Sampling Rate **/
#define LDAC_SMPLRATEBITS   (3)
/** Channel **/
#define LDAC_CHCONFIG2BITS  (2)
enum CHANNEL {
    MONO   = 0,
    STEREO = 1
};
/** Frame Length **/
#define LDAC_FRAMELEN2BITS  (9)
/** Frame Status **/
#define LDAC_FRAMESTATBITS  (2)

/** Band Info **/
#define LDAC_NBANDBITS      (4)
#define LDAC_BAND_OFFSET    (2)

/** Band **/
#define LDAC_MAXNBANDS        16

/* Flag */
#define LDAC_FLAGBITS       (1)
#define LDAC_TRUE           (1)
#define LDAC_FALSE          (0)

 /* Mode */
#define LDAC_MODE_0         (0)
#define LDAC_MODE_1         (1)
#define LDAC_MODE_2         (2)
#define LDAC_MODE_3         (3)

/** Gradient Data **/
#define LDAC_GRADMODEBITS      (2)
#define LDAC_GRADOSBITS        (5)
#define LDAC_GRADQU0BITS       (6)
#define LDAC_GRADQU1BITS       (5)
#define LDAC_NADJQUBITS        (5)

/** Scale Factor Data **/
#define LDAC_SFCMODEBITS       1
#define LDAC_IDSFBITS          5
#define LDAC_NSFCWTBL          8
#define LDAC_SFCBLENBITS       2
#define LDAC_MINSFCBLEN_0      3
#define LDAC_MINSFCBLEN_1      2
#define LDAC_MINSFCBLEN_2      2
#define LDAC_SFCWTBLBITS       3

#define LDAC_MINIDWL1          1
#define LDAC_MAXIDWL1         15
#define LDAC_MAXIDWL2         15

// Quantization Units
#define LDAC_MAXNQUS          34
#define LDAC_MAXGRADQU        50

/***************************************************************************************************
    Weighting Tables for Scale Factor Data
***************************************************************************************************/
static const uint8_t gaa_sfcwgt_ldac[LDAC_NSFCWTBL][LDAC_MAXNQUS] = {
{
     1,  0,  0,  1,  1,  1,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,
     3,  3,  3,  3,  3,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8,  8,
},
{
     0,  1,  1,  2,  3,  4,  4,  4,  4,  5,  6,  6,  6,  6,  6,  7,
     7,  7,  7,  7,  7,  7,  8,  8,  8,  9, 10, 10, 11, 11, 12, 12, 12, 12,
},
{
     0,  1,  1,  2,  3,  3,  3,  3,  3,  4,  4,  5,  5,  5,  5,  5,
     5,  5,  5,  5,  5,  5,  6,  6,  6,  7,  8,  9,  9, 10, 10, 11, 11, 11,
},
{
     0,  1,  3,  4,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  7,  7,
     7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  9,  9,  9, 10, 10, 10, 10,
},
{
     0,  1,  3,  4,  5,  5,  6,  7,  7,  8,  8,  9,  9, 10, 10, 10,
    10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13,
},
{
     1,  0,  1,  2,  2,  3,  3,  4,  4,  5,  6,  7,  7,  8,  8,  8,
     9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11,
},
{
     0,  0,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4,
     4,  4,  4,  4,  4,  4,  4,  5,  5,  6,  7,  7,  7,  8,  9,  9,  9,  9,
},
{
     0,  0,  1,  2,  3,  4,  4,  5,  5,  6,  7,  7,  8,  8,  8,  8,
     9,  9,  9,  9,  9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 12,
},
};

static const uint8_t ga_nqus_ldac[LDAC_MAXNBANDS+1] = {
     0,  4,  8, 10, 12, 14, 16, 18, 20, 22, 24, 25, 26, 28, 30, 32, 34,
};


int ldacdecInit( ldacdec_t *handle )
{
    InitHuffmanCodebooks();
    InitMdct();

    memset( &handle->frame.channels[0].mdct, 0, sizeof( Mdct ) ); 
    memset( &handle->frame.channels[1].mdct, 0, sizeof( Mdct ) );
    
    handle->frame.channels[0].frame = &handle->frame;
    handle->frame.channels[1].frame = &handle->frame;

    return 0;
}

static int decodeBand( frame_t *handle, BitReaderCxt *br )
{
    handle->nbrBands = ReadInt( br, LDAC_NBANDBITS ) + LDAC_BAND_OFFSET;
    LOG("nbrBands:        %d\n", handle->nbrBands );
    ReadInt( br, LDAC_FLAGBITS ); // unused

    handle->quantizationUnitCount = ga_nqus_ldac[handle->nbrBands];
    return 0;
}

static int decodeGradient( frame_t *handle, BitReaderCxt *br )
{
    handle->gradientMode = ReadInt( br, LDAC_GRADMODEBITS );
    if( handle->gradientMode == LDAC_MODE_0 )
    {
        handle->gradientStartUnit  = ReadInt( br, LDAC_GRADQU0BITS );
        handle->gradientEndUnit    = ReadInt( br, LDAC_GRADQU0BITS ) + 1;
        handle->gradientStartValue = ReadInt( br, LDAC_GRADOSBITS );
        handle->gradientEndValue   = ReadInt( br, LDAC_GRADOSBITS );
        LOG("gradient\n\tqu [%3d,%3d]\n\tvalue [%3d,%3d]\n", 
                handle->gradientStartUnit, handle->gradientEndUnit, 
                handle->gradientStartValue, handle->gradientEndValue );
    } else
    {
        handle->gradientStartUnit  = ReadInt( br, LDAC_GRADQU1BITS );
        handle->gradientEndUnit    = 26;
        handle->gradientStartValue = ReadInt( br, LDAC_GRADOSBITS );
        handle->gradientEndValue   = 31;
        LOG("gradient\n\tqu [%3d,%3d\n\tvalue [%3d,%3d]\n", 
                handle->gradientStartUnit, handle->gradientEndUnit, 
                handle->gradientStartValue, handle->gradientEndValue );
    }
    
    handle->gradientBoundary = ReadInt( br, LDAC_NADJQUBITS );
    return 0;
}

static void calculateGradient( frame_t *handle )
{
    int valueCount = handle->gradientEndValue - handle->gradientStartValue;
    int unitCount = handle->gradientEndUnit - handle->gradientStartUnit;
    
    memset( handle->gradient, 0, sizeof( handle->gradient ) );

    for( int i=0; i<handle->gradientEndUnit; ++i )
        handle->gradient[i] = -handle->gradientStartValue;
    for( int i=handle->gradientEndUnit; i<handle->quantizationUnitCount; ++i )
        handle->gradient[i] = -handle->gradientEndValue;

    if( unitCount > 0 )
    {
        const uint8_t *curve = gradientCurves[unitCount-1];
        if(  valueCount > 0 )
        {
            for( int i=handle->gradientStartUnit; i<handle->gradientEndUnit; ++i )
            {
                handle->gradient[i] -= ((curve[i-handle->gradientStartUnit] * (valueCount-1)) >> 8) + 1;
            }
        } else if( valueCount < 0 )
        {
            for( int i=handle->gradientStartUnit; i<handle->gradientEndUnit; ++i )
            {
                handle->gradient[i] -= ((curve[i-handle->gradientStartUnit] * (valueCount-1)) >> 8) + 1;
            }
        }
    }
    
    LOG_ARRAY_LEN( handle->gradient, "%3d, ", handle->quantizationUnitCount );
}

static void calculatePrecisionMask(channel_t* handle)
{
	memset(handle->precisionMask, 0, sizeof(handle->precisionMask));
	for (int i = 1; i < handle->frame->quantizationUnitCount; i++)
	{
		const int delta = handle->scaleFactors[i] - handle->scaleFactors[i - 1];
		if (delta > 1)
		{
			handle->precisionMask[i] += LDAC_MIN(delta - 1, 5);
		}
		else if (delta < -1)
		{
			handle->precisionMask[i - 1] += LDAC_MIN(delta * -1 - 1, 5);
		}
	}
}


static void calculatePrecisions( channel_t *handle )
{
    frame_t *frame = handle->frame;
    
    for( int i=0; i<frame->quantizationUnitCount; ++i )
    {
        switch( frame->gradientMode )
        {
            case LDAC_MODE_0:
            {
                int precision = handle->scaleFactors[i] + frame->gradient[i];
                precision = LDAC_MAX( precision, LDAC_MINIDWL1 );
                handle->precisions[i] = precision;
                break;
            }
            case LDAC_MODE_1:
            {
                int precision = handle->scaleFactors[i] + frame->gradient[i] + handle->precisionMask[i];
                if( precision > 0 )
                    precision /= 2;
                precision = LDAC_MAX( precision, LDAC_MINIDWL1 );
                handle->precisions[i] = precision;
                break;
            }
            case LDAC_MODE_2:
            {
                int precision = handle->scaleFactors[i] + frame->gradient[i] + handle->precisionMask[i];
                if( precision > 0 )
                    precision = ( precision * 3 ) / 8;
                precision = LDAC_MAX( precision, LDAC_MINIDWL1 );
                handle->precisions[i] = precision;
                break;
            }
            case LDAC_MODE_3:
            {
                int precision = handle->scaleFactors[i] + frame->gradient[i] + handle->precisionMask[i];
                if( precision > 0 )
                    precision /= 4;
                precision = LDAC_MAX( precision, LDAC_MINIDWL1 );
                handle->precisions[i] = precision;
                break;
            }
            default:
                assert(0);
                break;
        }
    }
    
    for( int i=0; i<frame->gradientBoundary; ++i )
    {
        handle->precisions[i]++;
    }

    for( int i=0; i<frame->quantizationUnitCount; ++i )
    {
        handle->precisionsFine[i] = 0;
        if( handle->precisions[i] > LDAC_MAXIDWL1 )
        {
            handle->precisionsFine[i] = handle->precisions[i] - LDAC_MAXIDWL1;
            handle->precisions[i] = LDAC_MAXIDWL1;
        }
    }

    LOG_ARRAY_LEN( handle->precisions, "%3d, ", frame->quantizationUnitCount );
    LOG_ARRAY_LEN( handle->precisionsFine, "%3d, ", frame->quantizationUnitCount );
}

static int decodeScaleFactor0( channel_t *handle, BitReaderCxt *br )
{
    LOG_FUNCTION();
    frame_t *frame = handle->frame;
    handle->scaleFactorBitlen = ReadInt( br, LDAC_SFCBLENBITS ) + LDAC_MINSFCBLEN_0;
    handle->scaleFactorOffset = ReadInt( br, LDAC_IDSFBITS );
    handle->scaleFactorWeight = ReadInt( br, LDAC_SFCWTBLBITS );

    LOG("scale factor bitlen = %d\n", handle->scaleFactorBitlen );
    LOG("scale factor offset = %d\n", handle->scaleFactorOffset );
    LOG("scale factor weight = %d\n", handle->scaleFactorWeight );

    const int mask = (1<<handle->scaleFactorBitlen)-1;
    const uint8_t *weightTable = gaa_sfcwgt_ldac[handle->scaleFactorWeight];
    const HuffmanCodebook* codebook = &HuffmanScaleFactorsUnsigned[handle->scaleFactorBitlen];
    handle->scaleFactors[0] = ReadInt( br, handle->scaleFactorBitlen );
//    printf("diff =\n    ");
    for( int i=1; i<frame->quantizationUnitCount; ++i )
    {
        int diff = ReadHuffmanValue( codebook, br, 1 );
       
//        printf("%2d, ", diff );
        handle->scaleFactors[i] = (handle->scaleFactors[i-1] + diff) & mask;
        handle->scaleFactors[i-1] += handle->scaleFactorOffset - weightTable[i-1]; // cancel weights
    }
    handle->scaleFactors[frame->quantizationUnitCount-1] += handle->scaleFactorOffset - weightTable[frame->quantizationUnitCount-1];

    LOG_ARRAY_LEN( handle->scaleFactors, "%2d, ", frame->quantizationUnitCount );
    return 0;
}

static int decodeScaleFactor1( channel_t *handle, BitReaderCxt *br )
{
    LOG_FUNCTION();
    frame_t *frame = handle->frame;
    handle->scaleFactorBitlen = ReadInt( br, LDAC_SFCBLENBITS ) + LDAC_MINSFCBLEN_1;
    
    if( handle->scaleFactorBitlen > 4 )
    {
        for( int i=0; i<frame->quantizationUnitCount; ++i )
        {
            handle->scaleFactors[i] = ReadInt( br, LDAC_IDSFBITS );
        }
    } else
    {
        handle->scaleFactorOffset = ReadInt( br, LDAC_IDSFBITS );
        handle->scaleFactorWeight = ReadInt( br, LDAC_SFCWTBLBITS );
        const uint8_t *weightTable = gaa_sfcwgt_ldac[handle->scaleFactorWeight];
        for( int i=0; i<frame->quantizationUnitCount; ++i )
        {
            handle->scaleFactors[i] = ReadInt( br, handle->scaleFactorBitlen ) - weightTable[i] + handle->scaleFactorOffset;
        }
    }
    return 0;
}

int decodeScaleFactor2( channel_t *handle, BitReaderCxt *br )
{
    LOG_FUNCTION();
    frame_t *frame = handle->frame;
    channel_t *other = &frame->channels[0];

    handle->scaleFactorBitlen = ReadInt( br, LDAC_SFCBLENBITS ) + LDAC_MINSFCBLEN_2;
    LOG("scale factor bitlen: %d\n", handle->scaleFactorBitlen );
    const HuffmanCodebook* codebook = &HuffmanScaleFactorsSigned[handle->scaleFactorBitlen];
//    printf("diff =\n");
    for( int i=0; i<frame->quantizationUnitCount; ++i )
    {
       int diff = ReadHuffmanValue( codebook, br, 1 );
//       printf("%2d, ", diff );
       handle->scaleFactors[i] = other->scaleFactors[i] + diff;
    }
    
    LOG_ARRAY_LEN( handle->scaleFactors, "%2d, ", frame->quantizationUnitCount );
    return 0;
}

int decodeScaleFactors( frame_t *handle, BitReaderCxt *br, int channelNbr )
{
    LOG_FUNCTION();
    channel_t *channel = &handle->channels[channelNbr];
    channel->scaleFactorMode = ReadInt( br, LDAC_SFCMODEBITS );
    LOG("scale factor mode = %d\n", channel->scaleFactorMode );
    if( channelNbr == 0 )
    {
        if( channel->scaleFactorMode == LDAC_MODE_0 )
            decodeScaleFactor0( channel, br );
        else
            decodeScaleFactor1( channel, br );
    } else
    {
        if( channel->scaleFactorMode == LDAC_MODE_0 )
            decodeScaleFactor0( channel, br );
        else
            decodeScaleFactor2( channel, br );
    }
    return 0;
}

enum {
    CHANNEL_1CH = 1,
    CHANNEL_2CH = 2,
};

static const char gaa_block_setting_ldac[4][4]=
{
    {CHANNEL_1CH, 1, MONO},
    {CHANNEL_2CH, 2, MONO, MONO},
    {CHANNEL_2CH, 1, STEREO},
    {0, 0, 0},
};

static void pcmFloatToShort( frame_t *handle, int16_t *pcmOut )
{
    int i=0;
    for(int smpl=0; smpl<handle->frameSamples; ++smpl )
    {
        for( int ch=0; ch<handle->channelCount; ++ch, ++i )
        {
            pcmOut[i] = Clamp16(Round(handle->channels[ch].pcm[smpl]));
        }
    }
}

static const int channelConfigIdToChannelCount[] = { 1, 2, 2 };

int ldacdecGetChannelCount( ldacdec_t *handle )
{
    return channelConfigIdToChannelCount[handle->frame.channelConfigId];
}

static const unsigned short sampleRateIdToSamplesPower[] = {
    7, 7, 8, 8
};

static const int sampleRateIdToFrequency[] = { 44100, 48000, 88200, 96000 };

int ldacdecGetSampleRate( ldacdec_t *handle )
{
    return sampleRateIdToFrequency[handle->frame.sampleRateId];
}

static int decodeFrame( frame_t *handle, BitReaderCxt *br )
{
    int syncWord = ReadInt( br, LDAC_SYNCWORDBITS );
    if( syncWord != LDAC_SYNCWORD )
        return -1;

    handle->sampleRateId = ReadInt( br, LDAC_SMPLRATEBITS );
    handle->channelConfigId = ReadInt( br, LDAC_CHCONFIG2BITS );
    handle->frameLength = ReadInt( br, LDAC_FRAMELEN2BITS ) + 1;
    handle->frameStatus = ReadInt( br, LDAC_FRAMESTATBITS );
    
    handle->channelCount = channelConfigIdToChannelCount[handle->channelConfigId];
    handle->frameSamplesPower = sampleRateIdToSamplesPower[handle->sampleRateId];
    handle->frameSamples = 1<<handle->frameSamplesPower;

    handle->channels[0].mdct.Bits = handle->frameSamplesPower;
    handle->channels[1].mdct.Bits = handle->frameSamplesPower;

    LOG("sampleRateId:    %d\n", handle->sampleRateId );
    LOG("   sample rate:  %d\n", sampleRateIdToFrequency[handle->sampleRateId] );
    LOG("   samplePower:  %d\n", handle->frameSamplesPower );
    LOG("channelConfigId: %d\n", handle->channelConfigId );
    LOG("frameLength:     %d\n", handle->frameLength );
    LOG("frameStatus:     %d\n", handle->frameStatus );

    return 0;
}

int ldacDecode( ldacdec_t *handle, uint8_t *stream, int16_t *pcm, int *bytesUsed )
{
    BitReaderCxt brObject;
    BitReaderCxt *br = &brObject;
    InitBitReaderCxt( br, stream );

    frame_t *frame = &handle->frame;
   
    int ret = decodeFrame( frame, br );
    if( ret < 0 )
        return -1;
   
    for( int block = 0; block<gaa_block_setting_ldac[frame->channelConfigId][1]; ++block )
    {
        decodeBand( frame, br );
        decodeGradient( frame, br );
        calculateGradient( frame );
        
        for( int i=0; i<frame->channelCount; ++i )
        {
            channel_t *channel = &frame->channels[i];
            decodeScaleFactors( frame, br, i );
            calculatePrecisionMask( channel ); 
            calculatePrecisions( channel );

            decodeSpectrum( channel, br );
            decodeSpectrumFine( channel, br );
            dequantizeSpectra( channel );
            scaleSpectrum( channel );

            RunImdct( &channel->mdct, channel->spectra, channel->pcm );
        }
        AlignPosition( br, 8 );

        pcmFloatToShort( frame, pcm );
    }
    AlignPosition( br, (frame->frameLength)*8 + 24 );

    if( bytesUsed != NULL )
        *bytesUsed = br->Position / 8;
    return 0;
}

// for packet loss concealment
static const int sa_null_data_size_ldac[2] = {
    11, 15,
};

static const uint8_t saa_null_data_ldac[2][15] = {
    {0x07, 0xa0, 0x16, 0x00, 0x20, 0xad, 0x51, 0x45, 0x14, 0x50, 0x49},
    {0x07, 0xa0, 0x0a, 0x00, 0x20, 0xad, 0x51, 0x41, 0x24, 0x93, 0x00, 0x28, 0xa0, 0x92, 0x49},
};

int ldacNullPacket( ldacdec_t *handle, uint8_t *output, int *bytesUsed )
{
    frame_t *frame = &handle->frame;
    uint8_t *ptr = output;

    for( int block = 0; block<gaa_block_setting_ldac[frame->channelConfigId][1]; ++block )
    {
        const int channelType = gaa_block_setting_ldac[frame->channelConfigId][2];
        const int size = sa_null_data_size_ldac[channelType];
        memcpy( ptr, saa_null_data_ldac[channelType], size );
        ptr += size;
    }

    *bytesUsed = frame->frameLength + 3;

    return 0;
}
