# Chapter 5: Configuring Model Input and Output

The reading and writing of model fields in MPAS is handled by user-configurable *streams*. A stream represents a fixed set of model fields, together with dimensions and attributes, that are all written or read together to or from the same file or set of files. Each MPAS model core may define its own set of default streams that it typically uses for reading initial conditions, for writing and reading restart fields, and for writing additional model history fields. Besides these default streams, users may define new streams to, e.g., write certain diagnostic fields at a higher temporal frequency than the usual model history fields.

Streams are defined in XML configuration files that are created at build time for each model core. The name of this XML file is simply `streams.` suffixed with the name of the core. For example, the streams for the *atmosphere* core are defined in a file named `streams.atmosphere`, and the streams for the *init_atmosphere* core are defined in a file named `streams.init_atmosphere`. An XML stream file may further reference other text files that contain lists of the model fields that are read or written in each of the streams defined in the XML stream file.

Changes to the XML stream configuration file will take effect the next time an MPAS core is run; there is no need to re-compile after making modifications to the XML files. It is therefore possible, e.g., to change the interval at which a stream is written, the template for the filenames associated with a stream, or the set of fields that are written to a stream, without the need to re-compile any code.

Two classes of streams exist in MPAS: *immutable* streams and *mutable* streams. Immutable streams are those for which the set of fields that belong to the stream may not be modified at model run-time; however, it is possible to modify the interval at which the stream is read or written, the filename template describing the files containing the stream on disk, and several other parameters of the stream. In contrast, all aspects of mutable streams, including the set of fields that belong to the stream, may be modified at run-time. The motivation for the creation of two stream classes is the idea that an MPAS core may not function correctly if certain fields are not read in upon model start-up or written to restart files, and it is therefore not reasonable for users to modify this set of required fields at run-time. An MPAS core developer may choose to implement such streams as immutable streams. Since fields may not be added to an immutable stream at run-time, new immutable streams may not be defined at run-time, and the only type of new stream that may be defined at run-time is the mutable stream type.

## 5.1 XML Stream Configuration Files

The XML stream configuration file for an MPAS core always has a parent XML element named `streams`, within which individual streams are defined:

```xml
<streams>

    ... one or more stream definitions ...

</streams>
```

Immutable streams are defined with the `immutable_stream` element, and mutable streams are defined with the `stream` element:

```xml
<immutable_stream name="initial_conditions"
                  type="input"
                  filename_template="init.nc"
                  input_interval="initial_only"
                  />

<stream name="history"
        type="output"
        filename_template="output.$Y-$M-$D_$h.$m.$s.nc"
        output_interval="6:00:00" >

    ... model fields belonging to this stream ...

</stream>
```

As shown in the example stream definitions above, both classes of stream have the following required attributes:

- **`name`** -- A unique name used to refer to the stream.

- **`type`** -- The type of stream, either `"input"`, `"output"`, `"input;output"`, or `"none"`. A stream may be both an input and an output stream (i.e., `"input;output"`) if, for example, it is read once at model start-up to provide initial conditions and thereafter written periodically to provide model checkpoints. A stream may be defined as neither input nor output (i.e., `"none"`) for the purposes of defining a set of fields for inclusion in other streams. Note that, for immutable streams, the type attribute may not be changed at run-time.

- **`filename_template`** -- The template for files that exist or will be created by the stream. The filename template may include any of the following variables, which are expanded based on the simulated time at which files are first created:
  - `$Y` -- Year
  - `$M` -- Month
  - `$D` -- Day of the month
  - `$d` -- Day of the year
  - `$h` -- Hour
  - `$m` -- Minute
  - `$s` -- Second

  A filename template may include either a relative or an absolute path, in which case MPAS will attempt to create any directories in the path that do not exist, subject to filesystem permissions.

- **`input_interval`** -- For streams that have type `"input"` or `"input;output"`, the interval, beginning at the model initial time, at which the stream will be read. Possible values include a time interval specification in the format `"YYYY-MM-DD_hh:mm:ss"`; the value `"initial_only"`, which specifies that the stream is read only once at the model initial time; or the value `"none"`, which specifies that the stream is not read during a model run.

- **`output_interval`** -- For streams that have type `"output"` or `"input;output"`, the interval, beginning at the model initial time, at which the stream will be written. Possible values include a time interval specification in the format `"YYYY-MM-DD_hh:mm:ss"`; the value `"initial_only"`, which specifies that the stream is written only once at the model initial time; or the value `"none"`, which specifies that the stream is not written during a model run.

Finally, the set of fields that belong to a mutable stream may be specified with any combination of the following elements. Note that, for immutable streams, no fields are specified at run-time in the XML configuration file.

- **`var`** -- Associates the specified variable with the stream. The variable may be any of those defined in an MPAS core's `Registry.xml` file, but may not include individual constituent arrays from a `var_array`.
- **`var_array`** -- Associates all constituent variables in a `var_array`, defined in an MPAS core's `Registry.xml` file, with the stream.
- **`var_struct`** -- Associates all variables in a `var_struct`, defined in an MPAS core's `Registry.xml` file, with the stream.
- **`stream`** -- Associates all explicitly associated fields in the specified stream with the stream; streams are not recursively included.
- **`file`** -- Associates all variables listed in the specified text file, with one field per line, with the stream.

## 5.2 Optional Stream Attributes

Besides the required attributes described in the preceding section, several additional, optional attributes may be added to the definition of a stream.

- **`filename_interval`** -- The interval between the timestamps used in the construction of the names of files associated with a stream. Possible values include a time interval specification in the format `"YYYY-MM-DD_hh:mm:ss"`; the value `"none"`, indicating that only one file containing all times is associated with the stream; the value `"input_interval"` that, for input type streams, indicates that each time to be read from the stream will come from a unique file; or the value `"output_interval"` that, for output type streams, indicates that each time to be written to the stream will go to a unique file whose name is based on the timestamp of the data being written. The default value is `"input_interval"` for input type streams and `"output_interval"` for output type streams. For streams of type `"input;output"`, the default filename interval is `"input_interval"` if the input interval is an interval (i.e., not `"initial_only"`), or `"output_interval"` otherwise. Refer to Section 5.3.1 for an example.

- **`reference_time`** -- A time that is an integral number of filename intervals from the timestamp of any file associated with the stream. The default value is the start time of the model simulation. Refer to Section 5.3.3 for an example.

- **`clobber_mode`** -- Specifies how a stream should handle attempts to write to a file that already exists. Possible values include:
  - `"overwrite"` -- The stream is allowed to overwrite records in existing files and to append new records to existing files; records not explicitly written to are left untouched.
  - `"truncate"` or `"replace_files"` -- The stream is allowed to overwrite existing files, which are first truncated to remove any existing records; this is equivalent to replacing any existing files with newly created files of the same name.
  - `"append"` -- The stream is only allowed to append new records to existing files; existing records may not be overwritten.
  - `"never_modify"` -- The stream is not allowed to modify existing files in any way.

  The default clobber mode for streams is `"never_modify"`. Refer to Section 5.3.2 for an example.

- **`precision`** -- The precision with which real-valued fields will be written or read in a stream. Possible values include `"single"` for 4-byte real values, `"double"` for 8-byte real values, or `"native"`, which specifies that real-valued fields will be written or read in whatever precision the MPAS core was compiled. The default value is `"native"`. Refer to Section 5.3.1 for an example.

- **`packages`** -- A list of packages attached to the stream. A stream will be active (i.e., read or written) only if at least one of the packages attached to it is active, or if no packages at all are attached. Package names are provided as a semi-colon-separated list. Note that packages may only be defined in an MPAS core's `Registry.xml` file at build time. By default, no packages are attached to a stream.

- **`io_type`** -- The underlying library and file format that will be used to read or write a stream. Possible values include:
  - `"pnetcdf"` -- Read/write the stream with classic large-file NetCDF files (CDF-2) using the ANL Parallel-NetCDF library.
  - `"pnetcdf,cdf5"` -- Read/write the stream with large-variable files (CDF-5) using the ANL Parallel-NetCDF library.
  - `"netcdf"` -- Read/write the stream with classic large-file NetCDF files (CDF-2) using the Unidata serial NetCDF library.
  - `"netcdf4"` -- Read/write the stream with HDF-5 files using the Unidata parallel NetCDF-4 library.

  Note that the PIO library must have been built with support for the selected `io_type`. By default, all input and output streams are read and written using the `"pnetcdf"` option.

## 5.3 Stream Definition Examples

This section provides several example streams that make use of the optional stream attributes described in Section 5.2. All examples are of output streams, since it is more likely that a user will need to write additional fields than to read additional fields, which a model would need to be aware of; however, the concepts illustrated here translate directly to input streams as well.

### 5.3.1 Example: A Single-Precision Output Stream with One Month of Data per File

In this example, the optional attribute specification `filename_interval="01-00_00:00:00"` is added to force a new output file to be created for the stream every month. Note that the general format for time interval specifications is `YYYY-MM-DD_hh:mm:ss`, where any leading terms can be omitted; in this case, the year part of the interval is omitted. To reduce the file size, the specification `precision="single"` is also added to force real-valued fields to be written as 4-byte floating-point values, rather than the default of 8 bytes.

```xml
<stream name="diagnostics"
        type="output"
        filename_template="diagnostics.$Y-$M.nc"
        filename_interval="01-00_00:00:00"
        precision="single"
        output_interval="6:00:00" >

    <var name="u10"/>
    <var name="v10"/>
    <var name="t2"/>
    <var name="q2"/>

</stream>
```

The only fields that will be written to this stream are the hypothetical 10-m diagnosed wind components, the 2-m temperature, and the 2-m specific humidity variables. Also, note that the filename template only includes the year and month from the model valid time; this can be problematic when the simulation starts in the middle of a month, and a solution for this problem is illustrated in the example of Section 5.3.3.

### 5.3.2 Example: Appending Records to Existing Output Files

By default, streams will never modify existing files whose filenames match the name of a file that would otherwise be written during the course of a simulation. However, when restarting a simulation that is expected to add more records to existing output files, it can be useful to instruct the MPAS I/O system to append these records, thereby modifying existing files. This may be accomplished with the `clobber_mode` attribute.

```xml
<stream name="diagnostics"
        type="output"
        filename_template="diagnostics.$Y-$M.nc"
        filename_interval="01-00_00:00:00"
        precision="single"
        clobber_mode="append"
        output_interval="6:00:00" >

    <var name="u10"/>
    <var name="v10"/>
    <var name="t2"/>
    <var name="q2"/>

</stream>
```

In general, if MPAS were to attempt to write a record at a time that already existed in an output file, a `clobber_mode` of `"append"` would not permit the write to take place, since this would modify existing data; in `"append"` mode, only new records may be added. However, due to a peculiarity in the implementation of the `"append"` clobber mode, it may be possible for an output file to contain duplicate times. This can happen when the first record that is appended to an existing file has a timestamp not matching any in the file, after which, any record that is written -- regardless of whether its timestamp matches one already in the file -- will be appended to the end of the file. This situation may arise, for example, when restarting a model simulation with a shorter `output_interval` than was used in the original model simulation with an MPAS core that does not write the first output time for restart runs.

### 5.3.3 Example: Referencing Filename Intervals to a Time Other Than the Start Time

The example stream of the previous sections creates a new file each month during the simulation, and the filenames contain only the year and month of the timestamp when the file was created. If a simulation begins at 00 UTC on the first day of a month, then each file in the diagnostic stream will contain only output times that fall within the month in the filename. However, if a simulation were to begin in the middle of a month -- for example, the month of June, 2014 -- the first diagnostics output file would have a filename of `diagnostics.2014-06.nc`, but rather than containing only output fields valid in June, it would contain all fields written between the middle of June and the middle of July, at which point one month of simulation would have elapsed, and a new output file, `diagnostics.2014-07.nc`, would be created.

In order to ensure that the file `diagnostics.2014-06.nc` contained only data from June 2014, the `reference_time` attribute may be added such that the day, hour, minute, and second in the date and time represent the first day of the month at 00 UTC. In this example, the year and month of the reference time are not important, since the purpose of the reference time here is to describe to MPAS that the monthly filename interval begins (i.e., is referenced to) the first day of the month.

```xml
<stream name="diagnostics"
        type="output"
        filename_template="diagnostics.$Y-$M.nc"
        filename_interval="01-00_00:00:00"
        reference_time="2014-01-01_00:00:00"
        precision="single"
        clobber_mode="append"
        output_interval="6:00:00" >

    <var name="u10"/>
    <var name="v10"/>
    <var name="t2"/>
    <var name="q2"/>

</stream>
```

In general, the components of a timestamp, `YYYY-MM-DD_hh:mm:ss`, that are less significant than (i.e., to the right of) those contained in a filename template are important in a `reference_time`. For example, with a `filename_template` that contained only the year, the month component of the `reference_time` would become important to identify the month of the year on which the yearly basis for filenames would begin.
