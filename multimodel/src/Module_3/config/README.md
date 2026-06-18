以下のファイル名と, 目的値の上下限値を設定してください.

- 共通の設定
  - output_prefix: 出力ファイルのプレフィックス
  - instance_file: インスタンスファイル
  - fringe_tree_file: フリンジツリーファイル

- input_data: 各物性の設定
  - property: 物性名
  - model: 予測モデル
  - target_value_lower_bound: 目的値の下限値
  - target_value_upper_bound: 目的値の上限値

  - LRの場合:
    Module2までで使用もしくは作成したファイルを指定してください.
    - LR_filename: LRファイル(モデルファイル)
    - desc_filename: 特徴ベクトルファイル
    - desc_norm_filename: 正規化特徴ベクトルファイル
    - fringe_filename: フリンジファイル
    - values_filename: 値ファイル
    - sdf_filename: sdfファイル

  - RFの場合:
    Module2までで使用もしくは作成したファイルを指定してください.
    - RF_filename: RFファイル(モデルファイル)
    - desc_filename: 特徴ベクトルファイル
    - desc_norm_selected_filename: 正規化特徴ベクトルファイル
    - fringe_filename: フリンジファイル
    - values_filename: 値ファイル
    - sdf_filename: sdfファイル

  - ANNの場合:
    Module2までで使用もしくは作成したファイルを指定してください.
    - biases_filename: バイアスファイル
    - weights_filename: 重みファイル
    - desc_filename: 特徴ベクトルファイル
    - desc_norm_selected_filename: 正規化特徴ベクトルファイル
    - fringe_filename: フリンジファイル
    - values_filename: 値ファイル
    - sdf_filename: sdfファイル

set the following file names and the upper and lower bounds of the target values.

- Common settings
  - output_prefix: Prefix of the output file
  - instance_file: Instance file
  - fringe_tree_file: Fringe tree file

- input_data: Settings for each property
  - property: Property name
  - model: Prediction model
  - target_value_lower_bound: Lower bound of the target value
  - target_value_upper_bound: Upper bound of the target value

  - For LR:
    Specify the files used or created by Module2.
    - LR_filename: LR file (model file)
    - desc_filename: Feature vector file
    - desc_norm_filename: Normalized feature vector file
    - fringe_filename: Fringe file
    - values_filename: Value file
    - sdf_filename: sdf file

  - For RF:
    Specify the files used or created by Module2.
    - RF_filename: RF file (model file)
    - desc_filename: Feature vector file
    - desc_norm_selected_filename: Normalized feature vector file
    - fringe_filename: Fringe file
    - values_filename: Value file
    - sdf_filename: sdf file

  - For ANN:
    Specify the files used or created by Module2.
    - biases_filename: Bias file
    - weights_filename: Weight file
    - desc_filename: Feature vector file
    - desc_norm_selected_filename: Normalized feature vector file
    - fringe_filename: Fringe file
    - values_filename: Value file
    - sdf_filename: sdf file
