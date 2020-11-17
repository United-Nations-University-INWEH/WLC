
// SurfaceWaterTool water detection algorithm (web application version)

/*
This tool utilizes percentiles to obtain a single mosaic image from the full stack of Landsat satellite images. 
From this, the (Modified) Normalized Difference Water Index [(M)NDWI] is calculated, over which a threshold is 
applied to classify pixels as water. To reduce the effect of vegetation and (hill)shadows, masks using NDVI and 
HAND are applied on the result. Optionally, defringe and/or cloud busting functions can be run over all images 
before the percentile reduction occurs.
*/
// ------------------------------------------------------------------------------------------------------- //
// Parameters
// ------------------------------------------------------------------------------------------------------- //

// default values for input parameters
var default_startdate    = '2019-04-01';
var default_enddate      = '2019-06-30';
var default_do_months    = false;
var default_pcnt_perm    = 40;
var default_pcnt_temp    = 8;
var default_water_thresh = 0.3;
var default_ndvi_thresh  = 0.6;
var default_hand_thresh  = 50;
var default_cloud_thresh = 80;
var default_defringe     = true;

var RedRiver = ee.FeatureCollection('users/hamidmehmoodUNU/Merged_2019_MB_Flood_Extent');

RedRiver = RedRiver.geometry();
Map.addLayer(RedRiver, {color: 'red'}, 'RedRiver');

// default value for export resolution
var default_export_res = 30;

// Landsat band names
var LC457_BANDS = ['B1',    'B1',   'B2',    'B3',  'B4',  'B5',    'B7'];
var LC8_BANDS   = ['B1',    'B2',   'B3',    'B4',  'B5',  'B6',    'B7'];
var STD_NAMES   = ['blue2', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2'];

// list of month names
var list_months = ['January', 'February', 'March', 'April', 'May', 'June', 'July',
                   'August', 'September', 'October', 'November', 'December'];

// visual parameters for percentile layers
var pcnt_viz_params = {min:0.06, max:0.5, gamma:1.5};

// colour rendering for water layers
var water_style = '\
<RasterSymbolizer>\
  <ColorMap extended="true" >\
    <ColorMapEntry color="#ffffff" quantity="0.0" label="-1"/>\
    <ColorMapEntry color="#9999ff" quantity="1.0" label="-1"/>\
    <ColorMapEntry color="#00008b" quantity="2.0" label="-1"/>\
  </ColorMap>\
</RasterSymbolizer>';

// default map view
var start_location = [105.0, 12.0];
var start_zoom     = 10;

// ------------------------------------------------------------------------------------------------------- //
// Functions
// ------------------------------------------------------------------------------------------------------- //

// update of tool with user-specified input parameters
function SurfaceWaterToolResults(inputParams) {
  
  // get input params
  var time_start   = ee.Date(inputParams.date_start);
  var time_end     = ee.Date(inputParams.date_end);
  var do_months    = inputParams.do_months;
  var defringe     = inputParams.defringe;
  var pcnt_perm    = parseFloat(inputParams.pcnt_perm);
  var pcnt_temp    = parseFloat(inputParams.pctn_temp);
  var water_thresh = parseFloat(inputParams.water_thresh);
  var ndvi_thresh  = parseFloat(inputParams.veg_thresh);
  var hand_thresh  = parseFloat(inputParams.hand_thresh);
  var cloud_thresh = parseFloat(inputParams.cloud_thresh);
  
  // filter Landsat image collections
  var images_l4 = filterImages(l4, geometry2, [time_start, time_end]);
  var images_l5 = filterImages(l5, geometry2, [time_start, time_end]);
  var images_l7 = filterImages(l7, geometry2, [time_start, time_end]);
  var images_l8 = filterImages(l8, geometry2, [time_start, time_end]);
  
  // cloud busting
  if (cloud_thresh > 0) {
   images_l4 = bustCloudsIC(images_l4, cloud_thresh);
   images_l5 = bustCloudsIC(images_l5, cloud_thresh);
   images_l7 = bustCloudsIC(images_l7, cloud_thresh);
   images_l8 = bustCloudsIC(images_l8, cloud_thresh);
  }
  
  // rename bands to common names that are the same for the different satellites
  images_l4 = images_l4.select(LC457_BANDS, STD_NAMES);
  images_l5 = images_l5.select(LC457_BANDS, STD_NAMES);
  images_l7 = images_l7.select(LC457_BANDS, STD_NAMES);
  images_l8 = images_l8.select(LC8_BANDS, STD_NAMES);
  
  // apply defringing of L5/L7
  if (defringe === true) {
    images_l5 = images_l5.map(defringeLandsat);
    images_l7 = images_l7.map(defringeLandsat);
	}
  
  // merge image collections
  var images = ee.ImageCollection(images_l4.merge(images_l5).merge(images_l7).merge(images_l8));
  // print(images);
  
  // run algorithm
  var resultsAlgorithm    = SurfaceWaterAlgorithm(images, pcnt_perm, pcnt_temp, water_thresh, ndvi_thresh, hand_thresh);
  
  var prcnt_img_permanent = resultsAlgorithm.pcnt_perm;
  var prcnt_img_temporary = resultsAlgorithm.pcnt_temp;
  var MNDWI_permanent     = resultsAlgorithm.mndwi_perm;
  var MNDWI_temporary     = resultsAlgorithm.mndwi_temp;
  var NDVI_permanent_pcnt = resultsAlgorithm.ndvi_perm;
  var NDVI_temporary_pcnt = resultsAlgorithm.ndvi_temp;
  var HAND_mask           = resultsAlgorithm.mask_hand;
  var NDVI_mask_permanent = resultsAlgorithm.mask_ndvi_perm;
  var NDVI_mask_temporary = resultsAlgorithm.mask_ndvi_temp;
  var water_complete      = resultsAlgorithm.water;
  
  // run algorithm for each month
  var water_monthly = ee.List([null]);
  var temp_images   = null;  // initialize variable outside of if-statement to get rid of code editor warnings
  var temp_results  = null;  // initialize variable outside of if-statement to get rid of code editor warnings
  if (do_months === true) {
    // get monthly images through mapping over calendarRange filter
    var getMonthlyWater = function(month_index) {
      temp_images  = images.filter(ee.Filter.calendarRange(month_index, month_index, 'month'));
      temp_results = SurfaceWaterAlgorithm(temp_images, pcnt_perm, pcnt_temp, water_thresh, ndvi_thresh, hand_thresh);
      return temp_results.water;
    };
    water_monthly = ee.List.sequence(1, 12).map(getMonthlyWater);
  }
  // print(water_monthly);
  
  return {
    'img_pcnt_perm': prcnt_img_permanent,
    'img_pcnt_temp': prcnt_img_temporary,
    'mndwi_perm': MNDWI_permanent,
    'mndwi_temp': MNDWI_temporary,
    'ndvi_perm': NDVI_permanent_pcnt,
    'ndvi_temp': NDVI_temporary_pcnt,
    'mask_hand': HAND_mask.updateMask(HAND_mask),
    'mask_ndvi_perm': NDVI_mask_permanent.updateMask(NDVI_mask_permanent),
    'mask_ndvi_temp': NDVI_mask_temporary.updateMask(NDVI_mask_temporary),
    'water': water_complete.updateMask(water_complete),
    'water_monthly': water_monthly
  };
}

// surface water detection algorithm
function SurfaceWaterAlgorithm(images, pcnt_perm, pcnt_temp, water_thresh, ndvi_thresh, hand_thresh) {
  
  // calculate percentile images
  var prcnt_img_permanent = images.reduce(ee.Reducer.percentile([pcnt_perm])).rename(STD_NAMES).clip(geometry2);
  var prcnt_img_temporary = images.reduce(ee.Reducer.percentile([pcnt_temp])).rename(STD_NAMES).clip(geometry2);

  // MNDWI
  var MNDWI_permanent = prcnt_img_permanent.normalizedDifference(['green', 'swir1']);
  var MNDWI_temporary = prcnt_img_temporary.normalizedDifference(['green', 'swir1']);

  // water
  var water_permanent = MNDWI_permanent.gt(water_thresh);
  var water_temporary = MNDWI_temporary.gt(water_thresh);
  
  // get HAND mask
  var HAND_mask = HAND.gt(hand_thresh);
  
  // get NDVI masks
  var NDVI_permanent_pcnt = prcnt_img_permanent.normalizedDifference(['nir', 'red']);
  var NDVI_temporary_pcnt = prcnt_img_temporary.normalizedDifference(['nir', 'red']);
  var NDVI_mask_permanent = NDVI_permanent_pcnt.gt(ndvi_thresh);
  var NDVI_mask_temporary = NDVI_temporary_pcnt.gt(ndvi_thresh);
  
  // combined NDVI and HAND masks
  var NDVI_and_HAND_mask_permanent = NDVI_mask_permanent.add(HAND_mask);
  var NDVI_and_HAND_mask_temporary = NDVI_mask_temporary.add(HAND_mask);
  
  // apply NDVI and HAND masks
  var water_permanent_NDVImasked = water_permanent.eq(1).and(NDVI_mask_permanent.eq(0));
  var water_permanent_HANDmasked = water_permanent.eq(1).and(HAND_mask.eq(0));
  var water_permanent_masked     = water_permanent.eq(1).and(NDVI_and_HAND_mask_permanent.eq(0));
  
  var water_temporary_NDVImasked = water_temporary.eq(1).and(NDVI_mask_temporary.eq(0));
  var water_temporary_HANDmasked = water_temporary.eq(1).and(HAND_mask.eq(0));
  var water_temporary_masked     = water_temporary.eq(1).and(NDVI_and_HAND_mask_temporary.eq(0));

  // single image with permanent and temporary water
  // var water_complete = water_permanent.add(water_temporary);
  var water_complete = water_permanent_masked.add(water_temporary_masked);
  
  return {
    'pcnt_perm': prcnt_img_permanent,
    'pcnt_temp': prcnt_img_temporary,
    'mndwi_perm': MNDWI_permanent,
    'mndwi_temp': MNDWI_temporary,
    'ndvi_perm': NDVI_permanent_pcnt,
    'ndvi_temp': NDVI_temporary_pcnt,
    'mask_hand': HAND_mask,
    'mask_ndvi_perm': NDVI_mask_permanent,//.updateMask(NDVI_mask_permanent),
    'mask_ndvi_temp': NDVI_mask_temporary,//.updateMask(NDVI_mask_temporary),
    'water': water_complete.updateMask(water_complete)
  };
}

// filter images
function filterImages(image_collection, bounds, dates) {
  return image_collection. filterBounds(bounds).filterDate(dates[0], dates[1].advance(1, 'day'));
}

// cloud busting
// (https://code.earthengine.google.com/63f075a9e212f6ed4770af44be18a4fe, Ian Housman and Carson Stam)
// adapted to be able to use within another function with ICs and/or cloud thresholds only defined there:
function bustCloudsIC(ic, cloud_thresh) {
  return ic.map(function(img) {
    var t = img;
    var cs = ee.Algorithms.Landsat.simpleCloudScore(img).select('cloud');
    var out = img.mask(img.mask().and(cs.lt(ee.Number(cloud_thresh))));
    return out.copyProperties(t);
	});
}

// defringe Landsat 5 and/or 7
// Defringe algorithm credits:
// Author:
// Bonnie Ruefenacht, PhD
// Senior Specialist
// RedCastle Resources, Inc.
// Working onsite at: 
// USDA Forest Service 
// Remote Sensing Applications Center (RSAC) 
// 2222 West 2300 South
// Salt Lake City, UT 84119
// Office: (801) 975-3828 
// Mobile: (801) 694-9215
// Email: bruefenacht@fs.fed.us
// RSAC FS Intranet website: http://fsweb.rsac.fs.fed.us/
// RSAC FS Internet website: http://www.fs.fed.us/eng/rsac/
// Purpose: Remove the fringes of landsat 5 and 7 scenes.
//
// Kernel for masking fringes found in L5 and L7 imagery
var k = ee.Kernel.fixed(41, 41, 
[[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]);
var fringeCountThreshold = 279;  // Define number of non null observations for pixel to not be classified as a fringe
function defringeLandsat(img) {
  var m   = img.mask().reduce(ee.Reducer.min());
  var sum = m.reduceNeighborhood(ee.Reducer.sum(), k, 'kernel');
  sum = sum.gte(fringeCountThreshold);
  img = img.mask(img.mask().and(sum));
  return img;
}

// get input parameters from UI
function getInputParameters() {
  var params_list   = input.widgets();
  var date_start    = params_list.get(1).widgets().get(1).getValue();
  var date_end      = params_list.get(2).widgets().get(1).getValue();
  var do_months     = params_list.get(3).widgets().get(1).getValue();
  var pcnt_perm     = params_list.get(4).widgets().get(1).getValue();
  var pctn_temp     = params_list.get(5).widgets().get(1).getValue();
  var water_thresh  = params_list.get(6).widgets().get(1).getValue();
  var veg_thresh    = params_list.get(7).widgets().get(1).getValue();
  var hand_thresh   = params_list.get(8).widgets().get(1).getValue();
  var cloud_thresh  = params_list.get(9).widgets().get(1).getValue();
  var defringe      = params_list.get(10).widgets().get(1).getValue();
  return {
    'date_start': date_start,
    'date_end': date_end,
    'do_months': do_months,
    'pcnt_perm': pcnt_perm,
    'pctn_temp': pctn_temp,
    'water_thresh': water_thresh,
    'veg_thresh': veg_thresh,
    'hand_thresh': hand_thresh,
    'cloud_thresh': cloud_thresh,
    'defringe': defringe,
  };
}

// reset input parameters to default values defined at top of script
function resetParamsToDefault() {
  var params_list = input.widgets();
  params_list.get(1).widgets().get(1).setValue(default_startdate);
  params_list.get(2).widgets().get(1).setValue(default_enddate);
  params_list.get(3).widgets().get(1).setValue(default_do_months);
  params_list.get(4).widgets().get(1).setValue(default_pcnt_perm);
  params_list.get(5).widgets().get(1).setValue(default_pcnt_temp);
  params_list.get(6).widgets().get(1).setValue(default_water_thresh);
  params_list.get(7).widgets().get(1).setValue(default_ndvi_thresh);
  params_list.get(8).widgets().get(1).setValue(default_hand_thresh);
  params_list.get(9).widgets().get(1).setValue(default_cloud_thresh);
  params_list.get(10).widgets().get(1).setValue(default_defringe);
}

// update all elements of the tool
function SurfaceWaterToolUpdate(inputParams, results) {
  // reset all map layers
  Map.layers().reset();
  // add background layer
  Map.addLayer(ee.Image().byte().paint(geometry2, 1000), {}, 'background', false);
  Map.layers().get(0).setOpacity(0.5);
  //add Red River Flood Extent - Hamid
  Map.addLayer(image, {palette: ['blue', 'green', 'red']}, 'Red River Flood Extent 2019');
  // add other layers
  Map.addLayer(results.img_pcnt_perm.select(['swir1', 'nir', 'green']), pcnt_viz_params, 'percentile (permanent)', false);
  Map.addLayer(results.img_pcnt_temp.select(['swir1', 'nir', 'green']), pcnt_viz_params, 'percentile (temporary)', false);
  Map.addLayer(results.mndwi_perm, {min:-1, max:1}, 'MNDWI (permanent)', false);
  Map.addLayer(results.mndwi_temp, {min:-1, max:1}, 'MNDWI (temporary)', false);
  
  // Map.addLayer(HAND.clip(geometry), {}, 'HAND', false);
  // Map.addLayer(results.ndvi_perm, {min:-1, max:1}, 'NDVI (permanent)', false);
  // Map.addLayer(results.ndvi_temp, {min:-1, max:1}, 'NDVI (temporary)', false);
  Map.addLayer(results.mask_hand.clip(geometry2), {palette:['000000']}, 'mask (HAND)', false);
  Map.addLayer(results.mask_ndvi_perm, {palette:['000000']}, 'mask (NDVI, permanent)', false);
  Map.addLayer(results.mask_ndvi_temp, {palette:['000000']}, 'mask (NDVI, temporary)', false);
  Map.addLayer(results.water.sldStyle(water_style), {}, 'water', false);
  if (inputParams.do_months === true) {
    for (var i=0; i<12; i++) {
      Map.addLayer(ee.Image(results.water_monthly.get(i)).sldStyle(water_style), {}, 'water (' + list_months[i] + ')', false);
    }
  }
  // show standard water layer if not using months (or last month if using months)
  Map.layers().get(Map.layers().length()-1).setShown(true);
  // add AoI and export bounds
  Map.addLayer(ee.Image().byte().paint(geometry2, 0, 2), {}, 'project bounds', false);
  Map.addLayer(ee.Image().byte().paint(geometry2, 0, 2), {}, 'current session bounds', false);
  Map.addLayer(ee.Image().byte().paint(export_geom, 0, 2), {}, 'export bounds', false);
}

// -------------------------------------------------------------------------------------------------------- //
// User Interface
// -------------------------------------------------------------------------------------------------------- //

// intro with descriptions
var intro = ui.Panel([
  // title
  ui.Label({
    value: 'Historical Flood Mapping Tool',
    style: {fontSize: '20px', fontWeight: 'bold'}
  }),
  // intro text
  ui.Label({
    value: "Change the input parameters below and press the 'Update map' button to visualize new results.\
            Data can also be exported to your Google Drive as a .tif file.",
    style: {fontSize: '12px', padding: '0px  0px 0px 0px'}
  }),
  // extra text
  ui.Label({
    value: "For questions please contact hamid.mehmood@unu.edu",
    style: {fontSize: '11px', padding: '0px  0px 0px 0px'}
  }),
  // line to separate intro from rest of UI panel
  // ui.Label('------------------------------------------------------------------------------------')
]);

// legend for water classes (uses padding to give the box height and width)
var legend = ui.Panel({widgets:[
  ui.Label({value: 'Legend (water):', style: {fontWeight: 'bold'}}),
  ui.Panel({widgets:[
    ui.Label({
      style: {
        backgroundColor: '#00008b',
        padding: '8px',
        margin: '0 0 4px 8px'
      }
    }),
    ui.Label({
      value: 'Permanent',
      style: {margin: '0 0 4px 6px'}
    })], layout: ui.Panel.Layout.Flow('horizontal')}),
  ui.Panel({widgets:[
    ui.Label({
      style: {
        backgroundColor: '#9999ff',
        padding: '8px',
        margin: '0 0 4px 8px'
      }
    }),
    ui.Label({
      value: 'Temporary',
      style: {margin: '0 0 4px 6px'}
    })], layout: ui.Panel.Layout.Flow('horizontal')})
  ], style:{margin: '0px 0px 10px 0px'}
});

var input = ui.Panel([
  ui.Label({value: 'Parameters:', style: {fontWeight: 'bold'}}),
  ui.Panel([ui.Label('Start date:'), ui.Textbox({value:default_startdate, style:{margin:'0px 28px', width:'100px'}})], 
           ui.Panel.Layout.flow('horizontal')),
  ui.Panel([ui.Label('End date:'), ui.Textbox({value:default_enddate, style:{margin:'0px 35px', width:'100px'}})], 
           ui.Panel.Layout.flow('horizontal')),
  ui.Panel([ui.Label('Monthly averages:'), ui.Checkbox({value:default_do_months, style:{margin:'7px 0px 0px 29px'}})],
           ui.Panel.Layout.flow('horizontal')),
  ui.Panel([ui.Label('Percentile permanent:'), ui.Textbox({value:default_pcnt_perm, style:{margin:'0px 7px', width:'50px'}})], 
           ui.Panel.Layout.flow('horizontal')),
  ui.Panel([ui.Label('Percentile temporary:'), ui.Textbox({value:default_pcnt_temp, style:{margin:'0px 10px', width:'50px'}})], 
           ui.Panel.Layout.flow('horizontal')),
  ui.Panel([ui.Label('Water index threshold:'), ui.Textbox({value:default_water_thresh, style:{margin:'0px 4px', width:'50px'}})], 
           ui.Panel.Layout.flow('horizontal')),
  ui.Panel([ui.Label('NDVI threshold:'), ui.Textbox({value:default_ndvi_thresh, style:{margin:'0px 46px', width:'50px'}})], 
           ui.Panel.Layout.flow('horizontal')),
  ui.Panel([ui.Label('HAND threshold:'), ui.Textbox({value:default_hand_thresh, style:{margin:'0px 38px', width:'50px'}})], 
           ui.Panel.Layout.flow('horizontal')),
  ui.Panel([ui.Label('Cloud busting:'), ui.Textbox({value:default_cloud_thresh, style:{margin:'0px 53px', width:'50px'}})], 
           ui.Panel.Layout.flow('horizontal')),
  ui.Panel([ui.Label('Defringe L5/L7:'), ui.Checkbox({value:default_defringe, style:{margin:'7px 0px 0px 45px'}})], 
           ui.Panel.Layout.flow('horizontal'))
]);

// buttons to reset parameters or update the map
var resetAndUpdateButtons = ui.Panel([
  ui.Button('Reset parameters', resetParamsToDefault),
  ui.Button('Update map', function() {
    var inputParams = getInputParameters();
    var results     = SurfaceWaterToolResults(inputParams);
    SurfaceWaterToolUpdate(inputParams, results);
  })
], ui.Panel.Layout.flow('horizontal'));

// Export functionality
var exportUI = ui.Panel([
  ui.Label({value: 'Export water layer(s):', style: {fontWeight: 'bold'}}),
  ui.Panel([ui.Label('Water layer:'), ui.Checkbox({value:true, style:{margin:'7px 0px 0px 68px'}})], 
           ui.Panel.Layout.flow('horizontal')),
  ui.Panel([ui.Label('Metadata (CSV):'), ui.Checkbox({value:true, style:{margin:'7px 0px 0px 39px'}})], 
           ui.Panel.Layout.flow('horizontal')),
  ui.Panel([ui.Label('Resolution [m]:'), ui.Textbox({value:default_export_res, style:{margin:'0px 49px', width:'50px'}})], 
           ui.Panel.Layout.flow('horizontal')),
  ui.Button({label:'Export', style:{margin:'5px 8px'}, onClick: function() {
            var inputParams     = getInputParameters();
            var results         = SurfaceWaterToolResults(inputParams);
            var export_scale    = exportUI.widgets().get(3).widgets().get(1).getValue();
            if (exportUI.widgets().get(1).widgets().get(1).getValue() === true) {
              Export.image.toDrive({
                image: results.water,
                description: 'SurfaceWaterTool_FullPeriod',
                // folder: ,
                // fileNamePrefix: ,
                // dimensions: ,
                region: export_geom,
                scale: export_scale,
                // crs: ,
                // crsTransform: ,
                // maxPixels: 
              });
            }
            if (exportUI.widgets().get(2).widgets().get(1).getValue() === true) {
              Export.table.toDrive({
                collection: ee.FeatureCollection(ee.Feature(null, inputParams)),
                description: 'SurfaceWaterTool_Metadata',
                // folder: ,
                // fileNamePrefix: ,
                fileFormat: 'CSV'
              });
            }
          }}),
  ui.Label({value: 'Pressing the Export button queues up one or more exports in the Tasks tab\
    of the code editor, which need to be started manually.', style: {fontSize: '10px'}})
]);

// Panel combining all UI elements for the ui.root
var panel = ui.Panel({
  widgets: [intro, legend, input, resetAndUpdateButtons, exportUI],
  layout: ui.Panel.Layout.flow('vertical'),
  style: {
    position: 'top-left',
    width: '380px'
  }
});

// add the panel to the ui.root
ui.root.insert(0, panel);

// -------------------------------------------------------------------------------------------------------- //
// Initialization
// -------------------------------------------------------------------------------------------------------- //

// check if user imported geometry
if (geometry) {
  // if so, just center map on geometry
  Map.centerObject(geometry);
} else {
  // if not, use AoI instead
  var geometry = AOI.geometry();
  // center map on specified starting location
  Map.centerObject(ee.Geometry.Point(start_location), start_zoom);
}
var export_geom = geometry.bounds();

// assign large (positive!) HAND value to value found in strange horizontal lines (-99999),
// so it is masked unless user specifies a very large threshold
// (commented out for now)
// var HAND = HAND.where(HAND.lt(0), 1000);

// calculate initial results based on default parameters
var inputParams = getInputParameters();
var results     = SurfaceWaterToolResults(inputParams);

// show initial results
SurfaceWaterToolUpdate(inputParams, results);
